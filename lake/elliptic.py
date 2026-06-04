"""
Singular elliptic reconstruction for the invisc grad_perp psi.Singular elliptic reconstruction for the inviscid lake equations.

For numerical solution, we solve the equivalent symmetric positive definite
system:

    A_sing psi = rhs,

where:

    A_sing = -div((1/b) grad),
    rhs    = -b omega.

The coefficient 1/b is singular near the boundary because b -> 0.
This module preserves the coefficient 1/b in divergence form.

Boundary condition:

    psi = 0

on the discrete mask boundary.

This is a masked Cartesian prototype, not a boundary-fitted FEM/FVM
implementation.
"""

from dataclasses import dataclass
from typing import Dict, Any, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .config import LakeConfig
from .grid import LakeGrid
from .operators import (
    grad_perp_scalar,
    curl_vector,
    div_bu,
    masked_l2_norm,
    weighted_l2b_velocity_norm,
    ordinary_l2_vector_norm,
)


@dataclass
class SingularEllipticData:
    """
    Container for the singular elliptic reconstruction operator.

    Attributes
    ----------
    A_sing:
        Sparse matrix representing:

            A_sing = -div((1/b) grad)

        on the interior unknowns.

    idx:
        Integer index map from grid nodes to unknown indices.

    num_unknowns:
        Number of elliptic unknowns.

    symmetry_defect:
        Relative symmetry defect ||A-A.T||/||A||.

    nnz:
        Number of nonzero entries.
    """

    A_sing: sp.csr_matrix
    idx: np.ndarray
    num_unknowns: int
    symmetry_defect: float
    nnz: int


@dataclass
class ReconstructionResult:
    """
    Result of stream-function reconstruction from vorticity.

    Attributes
    ----------
    psi:
        Reconstructed stream function.

    u:
        Reconstructed velocity.

    info:
        Dictionary of solver diagnostics.
    """

    psi: np.ndarray
    u: np.ndarray
    info: Dict[str, Any]


# =============================================================================
# Sparse CG helper
# =============================================================================

def cg_solve_compat(
    A: sp.spmatrix,
    rhs: np.ndarray,
    rtol: float,
    atol: float,
    maxiter: int,
) -> Tuple[np.ndarray, int, int]:
    """
    Solve a sparse linear system with CG, supporting both recent and older SciPy.

    Recent SciPy versions use:

        cg(A, b, rtol=..., atol=...)

    Older versions use:

        cg(A, b, tol=...)

    Parameters
    ----------
    A:
        Sparse SPD matrix.

    rhs:
        Right-hand side.

    rtol, atol:
        Solver tolerances.

    maxiter:
        Maximum number of iterations.

    Returns
    -------
    solution:
        CG solution.

    info:
        CG info flag.

    iterations:
        Number of CG iterations counted by callback.
    """
    counter = {"count": 0}

    def callback(_xk):
        counter["count"] += 1

    try:
        sol, info = spla.cg(
            A,
            rhs,
            rtol=rtol,
            atol=atol,
            maxiter=maxiter,
            callback=callback,
        )
    except TypeError:
        sol, info = spla.cg(
            A,
            rhs,
            tol=rtol,
            maxiter=maxiter,
            callback=callback,
        )

    return sol, info, counter["count"]


# =============================================================================
# Face coefficient helper
# =============================================================================

def harmonic_face_average(
    a_left: float,
    a_right: float,
    eps: float = 1.0e-14,
) -> float:
    """
    Harmonic average for positive diffusion coefficients.

    For a divergence-form operator:

        -div(a grad psi),

    harmonic averaging is often more robust than arithmetic averaging when
    the coefficient a has strong spatial contrast.

    In the inviscid lake reconstruction:

        a = 1/b,

    and a becomes large near the degenerate boundary where b -> 0.
    Harmonic averaging avoids excessive overestimation of the face coefficient
    when neighboring values differ strongly.

    Parameters
    ----------
    a_left, a_right:
        Positive coefficients on neighboring nodes.

    eps:
        Safety value to avoid division by zero.

    Returns
    -------
    a_face:
        Harmonic face coefficient.
    """
    denom = a_left + a_right

    if denom <= eps:
        return 0.0

    return 2.0 * a_left * a_right / max(denom, eps)


# =============================================================================
# Assembly of A_sing = -div((1/b) grad)
# =============================================================================

def build_unknown_index(
    grid: LakeGrid,
) -> np.ndarray:
    """
    Build the unknown index map for interior elliptic unknowns.

    Unknowns are active nodes excluding the topological boundary band.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    Returns
    -------
    idx:
        Integer index map. Values are -1 outside unknowns.
    """
    idx = -np.ones_like(grid.mask, dtype=int)
    idx[grid.interior_unknowns] = np.arange(np.sum(grid.interior_unknowns))

    return idx


def assemble_singular_matrix(
    grid: LakeGrid,
) -> SingularEllipticData:
    """
    Assemble the singular elliptic matrix:

        A_sing = -div((1/b) grad)

    on the interior unknowns.

    The coefficient is:

        a = 1/b.

    Face coefficients between two active mask nodes are computed using
    harmonic averaging:

        a_face = 2 a_i a_j / (a_i + a_j).

    Harmonic averaging is more robust for high-contrast coefficients such as
    a=1/b near the degenerate boundary.

    Dirichlet boundary condition:

        psi = 0

    is imposed on the discrete mask boundary by adding only diagonal
    contributions when a neighbor is outside the set of unknowns.

    Boundary-face contributions are treated as first-order masked-grid flux
    approximations. This is consistent with the prototype masked Cartesian
    setting but is not a boundary-fitted FEM/FVM treatment.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    Returns
    -------
    data:
        SingularEllipticData object.
    """
    idx = build_unknown_index(grid)
    num_unknowns = int(np.sum(grid.interior_unknowns))

    a_coeff = np.zeros_like(grid.b_safe)
    a_coeff[grid.mask] = 1.0 / grid.b_safe[grid.mask]

    rows = []
    cols = []
    vals = []

    ny, nx = grid.mask.shape
    h2 = grid.h**2

    for i, j in np.argwhere(grid.interior_unknowns):
        p = idx[i, j]
        diag = 0.0

        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ii = i + di
            jj = j + dj

            if ii < 0 or ii >= ny or jj < 0 or jj >= nx:
                # Outside array: Dirichlet psi = 0.
                # First-order boundary flux approximation.
                a_face = a_coeff[i, j]
                diag += a_face / h2
                continue

            if grid.mask[ii, jj]:
                a_face = harmonic_face_average(
                    a_coeff[i, j],
                    a_coeff[ii, jj],
                )

                if grid.interior_unknowns[ii, jj]:
                    q = idx[ii, jj]
                    diag += a_face / h2
                    rows.append(p)
                    cols.append(q)
                    vals.append(-a_face / h2)
                else:
                    # Neighbor is boundary-band node: psi = 0.
                    # First-order boundary flux approximation.
                    diag += a_face / h2
            else:
                # Neighbor outside mask: psi = 0.
                # First-order boundary flux approximation using the interior
                # coefficient as the boundary-face coefficient.
                a_face = a_coeff[i, j]
                diag += a_face / h2

        rows.append(p)
        cols.append(p)
        vals.append(diag)

    A = sp.csr_matrix(
        (vals, (rows, cols)),
        shape=(num_unknowns, num_unknowns),
    )

    A_norm = float(spla.norm(A))
    symmetry_defect = float(spla.norm(A - A.T)) / max(A_norm, 1.0e-14)

    return SingularEllipticData(
        A_sing=A,
        idx=idx,
        num_unknowns=num_unknowns,
        symmetry_defect=symmetry_defect,
        nnz=A.nnz,
    )


def print_singular_matrix_summary(
    data: SingularEllipticData,
) -> None:
    """
    Print a summary of the singular elliptic matrix.

    Parameters
    ----------
    data:
        SingularEllipticData instance.
    """
    print("=" * 70)
    print("Singular elliptic matrix summary")
    print("=" * 70)
    print("Operator: A_sing = -div((1/b) grad)")
    print("Face averaging: harmonic average for interior mask faces")
    print("Boundary treatment: first-order masked-grid Dirichlet flux")
    print(f"shape              = {data.A_sing.shape}")
    print(f"number of unknowns = {data.num_unknowns}")
    print(f"nnz                = {data.nnz}")
    print(f"symmetry defect    = {data.symmetry_defect:.8e}")
    print("=" * 70)


# =============================================================================
# Reconstruction from vorticity
# =============================================================================

def reconstruct_from_vorticity(
    omega: np.ndarray,
    grid: LakeGrid,
    elliptic_data: SingularEllipticData,
    cfg: LakeConfig,
) -> ReconstructionResult:
    """
    Reconstruct stream function and velocity from vorticity.

    Continuous equation:

        div((1/b) grad psi) = b omega.

    SPD form solved numerically:

        A_sing psi = -b omega,

    where:

        A_sing = -div((1/b) grad).

    Velocity reconstruction:

        u = (1/b) grad_perp psi.

    Parameters
    ----------
    omega:
        Vorticity field.

    grid:
        LakeGrid instance.

    elliptic_data:
        SingularEllipticData containing A_sing and unknown index map.

    cfg:
        LakeConfig instance.

    Returns
    -------
    result:
        ReconstructionResult object.
    """
    A = elliptic_data.A_sing

    rhs_grid = -grid.b * omega
    rhs_vec = rhs_grid[grid.interior_unknowns]

    psi_vec, cg_info, cg_iterations = cg_solve_compat(
        A=A,
        rhs=rhs_vec,
        rtol=cfg.cg_rtol,
        atol=cfg.cg_atol,
        maxiter=cfg.cg_maxiter,
    )

    psi = np.zeros_like(omega)
    psi[grid.interior_unknowns] = psi_vec
    psi[grid.boundary_band_topological] = 0.0
    psi[~grid.mask] = 0.0

    u = grad_perp_scalar(
        psi,
        grid.dx,
        grid.dy,
        grid.mask,
    ) / grid.b_safe[..., None]

    u[~grid.mask] = 0.0

    residual_vec = A @ psi_vec - rhs_vec

    elliptic_residual = float(np.linalg.norm(residual_vec))
    rhs_norm = float(np.linalg.norm(rhs_vec))
    relative_elliptic_residual = elliptic_residual / max(rhs_norm, 1.0e-14)

    info = {
        "cg_info": cg_info,
        "cg_iterations": cg_iterations,
        "elliptic_residual": elliptic_residual,
        "relative_elliptic_residual": relative_elliptic_residual,
        "rhs_norm": rhs_norm,
    }

    return ReconstructionResult(
        psi=psi,
        u=u,
        info=info,
    )


# =============================================================================
# Reconstruction diagnostics
# =============================================================================

def reconstruction_diagnostics(
    psi_reference: np.ndarray,
    u_reference: np.ndarray,
    omega: np.ndarray,
    reconstruction: ReconstructionResult,
    grid: LakeGrid,
) -> Dict[str, float]:
    """
    Compute diagnostics for elliptic reconstruction.

    Parameters
    ----------
    psi_reference:
        Reference stream function, typically psi0.

    u_reference:
        Reference velocity, typically u0.

    omega:
        Input vorticity.

    reconstruction:
        ReconstructionResult.

    grid:
        LakeGrid instance.

    Returns
    -------
    diagnostics:
        Dictionary of reconstruction diagnostics.
    """
    psi_rec = reconstruction.psi
    u_rec = reconstruction.u

    psi_error = psi_rec - psi_reference
    u_error = u_rec - u_reference

    psi_ref_l2 = masked_l2_norm(
        psi_reference,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    psi_error_l2 = masked_l2_norm(
        psi_error,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    relative_psi_error = psi_error_l2 / max(psi_ref_l2, 1.0e-14)

    u_ref_l2b = weighted_l2b_velocity_norm(
        u_reference,
        grid.b,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    u_error_l2b = weighted_l2b_velocity_norm(
        u_error,
        grid.b,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    relative_velocity_error = u_error_l2b / max(u_ref_l2b, 1.0e-14)

    curl_u_rec = curl_vector(
        u_rec,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    curl_over_b_rec = curl_u_rec / grid.b_safe

    vorticity_residual = curl_over_b_rec - omega

    vorticity_residual_l2 = masked_l2_norm(
        vorticity_residual,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    omega_l2 = masked_l2_norm(
        omega,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    relative_vorticity_residual = vorticity_residual_l2 / max(
        omega_l2,
        1.0e-14,
    )

    div_bu_rec = div_bu(
        u_rec,
        grid.b,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    div_bu_rec_l2 = masked_l2_norm(
        div_bu_rec,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    bu_rec_l2 = ordinary_l2_vector_norm(
        grid.b[..., None] * u_rec,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    relative_div_bu_rec = div_bu_rec_l2 / max(
        bu_rec_l2 / grid.h,
        1.0e-14,
    )

    diagnostics = {
        "psi_error_l2": psi_error_l2,
        "relative_psi_error": relative_psi_error,
        "u_error_l2b": u_error_l2b,
        "relative_velocity_error": relative_velocity_error,
        "vorticity_residual_l2": vorticity_residual_l2,
        "relative_vorticity_residual": relative_vorticity_residual,
        "div_bu_rec_l2": div_bu_rec_l2,
        "relative_div_bu_rec": relative_div_bu_rec,
        "elliptic_residual": reconstruction.info["elliptic_residual"],
        "relative_elliptic_residual": reconstruction.info[
            "relative_elliptic_residual"
        ],
        "cg_info": reconstruction.info["cg_info"],
        "cg_iterations": reconstruction.info["cg_iterations"],
    }

    return diagnostics


def print_reconstruction_diagnostics(
    diagnostics: Dict[str, float],
) -> None:
    """
    Print diagnostics for elliptic reconstruction.

    Parameters
    ----------
    diagnostics:
        Dictionary returned by reconstruction_diagnostics.
    """
    print("=" * 70)
    print("Singular elliptic reconstruction diagnostics")
    print("=" * 70)

    print(f"||psi_rec - psi_ref||_L2              = {diagnostics['psi_error_l2']:.8e}")
    print(f"relative psi error                    = {diagnostics['relative_psi_error']:.8e}")
    print(f"||u_rec - u_ref||_L2_b                = {diagnostics['u_error_l2b']:.8e}")
    print(f"relative velocity error               = {diagnostics['relative_velocity_error']:.8e}")
    print(f"||curl(u_rec)/b - omega||_L2          = {diagnostics['vorticity_residual_l2']:.8e}")
    print(f"relative vorticity residual           = {diagnostics['relative_vorticity_residual']:.8e}")
    print(f"||div(bu_rec)||_L2                    = {diagnostics['div_bu_rec_l2']:.8e}")
    print(f"relative div(bu_rec)                  = {diagnostics['relative_div_bu_rec']:.8e}")
    print(f"elliptic residual                     = {diagnostics['elliptic_residual']:.8e}")
    print(f"relative elliptic residual            = {diagnostics['relative_elliptic_residual']:.8e}")
    print(f"CG info                               = {diagnostics['cg_info']}")
    print(f"CG iterations                         = {diagnostics['cg_iterations']}")

    print("=" * 70)


def reconstruction_sign_warning(
    psi_reference: np.ndarray,
    psi_reconstructed: np.ndarray,
    grid: LakeGrid,
) -> bool:
    """
    Detect a possible sign convention issue.

    If psi_reconstructed is closer to -psi_reference than to psi_reference,
    a sign warning is returned.

    Parameters
    ----------
    psi_reference:
        Reference stream function.

    psi_reconstructed:
        Reconstructed stream function.

    grid:
        LakeGrid instance.

    Returns
    -------
    warning:
        True if a possible sign issue is detected.
    """
    err_plus = masked_l2_norm(
        psi_reconstructed + psi_reference,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    err_minus = masked_l2_norm(
        psi_reconstructed - psi_reference,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    return err_plus < err_minus



