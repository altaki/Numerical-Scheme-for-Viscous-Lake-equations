"""
Degenerate weighted projection for the viscous lake equations.

This module implements the treatment:This module implements the weighted projection used in the viscous velocity

    phi = 0

on the discrete mask boundary.

Important prototype limitation
------------------------------
The sparse operator A_proj is assembled in divergence form, but the velocity
correction

    u_projected = u_star - grad(phi)

uses the lightweight finite-difference gradient grad_scalar from operators.py.
Therefore, the sparse divergence-form solve and the correction gradient are
not an exactly compatible algebraic gradient-divergence pair.

Consequently, exact discrete cancellation of div(b u_projected) is not
guaranteed, even when epsilon_proj = 0. This is why projection defects are
computed and reported.

This is a masked Cartesian prototype. It does not implement the full Navier
slip-with-friction boundary condition.
"""

from dataclasses import dataclass
from typing import Dict, Any, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .config import LakeConfig
from .grid import LakeGrid
from .operators import (
    grad_scalar,
    div_bu,
    masked_l2_norm,
    ordinary_l2_vector_norm,
)


@dataclass
class ProjectionData:
    """
    Container for the degenerate projection operator.

    Attributes
    ----------
    A_proj:
        Sparse matrix representing:

            A_proj = -div(b_proj grad)

        on projection unknowns.

    idx_proj:
        Integer index map for projection unknowns.

    projection_unknowns:
        Boolean mask of projection unknowns.

    b_proj:
        Projection coefficient.

    epsilon_proj:
        Stabilization parameter.

    num_unknowns:
        Number of projection unknowns.

    symmetry_defect:
        Relative symmetry defect ||A-A.T||/||A||.

    nnz:
        Number of nonzero matrix entries.
    """

    A_proj: sp.csr_matrix
    idx_proj: np.ndarray
    projection_unknowns: np.ndarray
    b_proj: np.ndarray
    epsilon_proj: float
    num_unknowns: int
    symmetry_defect: float
    nnz: int


@dataclass
class ProjectionResult:
    """
    Result of weighted projection.

    Attributes
    ----------
    u_projected:
        Projected velocity.

    phi:
        Projection potential.

    info:
        Dictionary of diagnostics.
    """

    u_projected: np.ndarray
    phi: np.ndarray
    info: Dict[str, Any]


# =============================================================================
# CG compatibility helper
# =============================================================================

def cg_solve_compat(
    A: sp.spmatrix,
    rhs: np.ndarray,
    rtol: float,
    atol: float,
    maxiter: int,
) -> Tuple[np.ndarray, int, int]:
    """
    Solve A x = rhs with conjugate gradients while supporting multiple SciPy
    versions.

    Recent SciPy versions use:

        cg(A, b, rtol=..., atol=...)

    Older SciPy versions use:

        cg(A, b, tol=...)

    Parameters
    ----------
    A:
        Sparse SPD matrix.

    rhs:
        Right-hand side vector.

    rtol, atol:
        Relative and absolute CG tolerances.

    maxiter:
        Maximum number of CG iterations.

    Returns
    -------
    sol:
        Solution vector.

    info:
        CG info flag.

    iterations:
        Number of iterations counted by callback.
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
# Projection unknowns and matrix assembly
# =============================================================================

def build_projection_unknowns(
    grid: LakeGrid,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build the projection unknown mask and index map.

    We use the same simple boundary treatment as for the singular elliptic
    problem:

        phi = 0 on the topological mask boundary.

    Unknowns are:

        projection_unknowns = mask & (~boundary_band_topological).

    Parameters
    ----------
    grid:
        LakeGrid instance.

    Returns
    -------
    projection_unknowns:
        Boolean mask.

    idx_proj:
        Integer index map with -1 outside unknowns.
    """
    projection_unknowns = grid.mask & (~grid.boundary_band_topological)

    idx_proj = -np.ones_like(grid.mask, dtype=int)
    idx_proj[projection_unknowns] = np.arange(np.sum(projection_unknowns))

    return projection_unknowns, idx_proj


def build_projection_coefficient(
    grid: LakeGrid,
    epsilon_proj: float,
) -> np.ndarray:
    """
    Build the projection coefficient b_proj.

    If epsilon_proj > 0:

        b_proj = b + epsilon_proj

    inside the mask.

    If epsilon_proj = 0:

        b_proj = b.

    Outside the mask, b_proj = 0.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    epsilon_proj:
        Projection stabilization parameter.

    Returns
    -------
    b_proj:
        Projection coefficient.
    """
    b_proj = np.zeros_like(grid.b)

    if epsilon_proj > 0.0:
        b_proj[grid.mask] = grid.b[grid.mask] + epsilon_proj
    else:
        b_proj[grid.mask] = grid.b[grid.mask]

    b_proj[~grid.mask] = 0.0

    return b_proj


def assemble_degenerate_projection_matrix(
    grid: LakeGrid,
    b_proj: np.ndarray,
    projection_unknowns: np.ndarray,
    idx_proj: np.ndarray,
    epsilon_proj: float,
) -> ProjectionData:
    """
    Assemble the degenerate projection matrix:

        A_proj = -div(b_proj grad)

    on projection_unknowns.

    Five-point divergence-form stencil with arithmetic face averaging of
    b_proj.

    Since b_proj = b + epsilon_proj is positive inside the mask when
    regularization is used, arithmetic averaging is acceptable for this
    prototype. More structure-preserving alternatives could be tested in a
    research-grade implementation.

    Boundary treatment:

        phi = 0

    on the discrete mask boundary.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    b_proj:
        Projection coefficient.

    projection_unknowns:
        Boolean mask of unknowns.

    idx_proj:
        Integer index map.

    epsilon_proj:
        Stabilization parameter.

    Returns
    -------
    data:
        ProjectionData object.
    """
    rows = []
    cols = []
    vals = []

    ny, nx = grid.mask.shape
    h2 = grid.h**2

    for i, j in np.argwhere(projection_unknowns):
        p = idx_proj[i, j]
        diag = 0.0

        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            ii = i + di
            jj = j + dj

            if ii < 0 or ii >= ny or jj < 0 or jj >= nx:
                # Outside array: Dirichlet phi = 0.
                # First-order boundary flux approximation.
                b_face = b_proj[i, j]
                diag += b_face / h2
                continue

            if grid.mask[ii, jj]:
                b_face = 0.5 * (b_proj[i, j] + b_proj[ii, jj])

                if projection_unknowns[ii, jj]:
                    q = idx_proj[ii, jj]
                    diag += b_face / h2
                    rows.append(p)
                    cols.append(q)
                    vals.append(-b_face / h2)
                else:
                    # Neighbor is on boundary band: phi = 0.
                    # First-order boundary flux approximation.
                    diag += b_face / h2
            else:
                # Neighbor outside mask: phi = 0.
                # First-order boundary flux approximation.
                b_face = b_proj[i, j]
                diag += b_face / h2

        rows.append(p)
        cols.append(p)
        vals.append(diag)

    num_unknowns = int(np.sum(projection_unknowns))

    A_proj = sp.csr_matrix(
        (vals, (rows, cols)),
        shape=(num_unknowns, num_unknowns),
    )

    A_norm = float(spla.norm(A_proj))
    symmetry_defect = float(spla.norm(A_proj - A_proj.T)) / max(A_norm, 1.0e-14)

    return ProjectionData(
        A_proj=A_proj,
        idx_proj=idx_proj,
        projection_unknowns=projection_unknowns,
        b_proj=b_proj,
        epsilon_proj=float(epsilon_proj),
        num_unknowns=num_unknowns,
        symmetry_defect=symmetry_defect,
        nnz=A_proj.nnz,
    )


def build_projection_data(
    grid: LakeGrid,
    cfg: LakeConfig,
    epsilon_proj: float | None = None,
) -> ProjectionData:
    """
    Build all projection data.

    If epsilon_proj is None, the default is:

        epsilon_proj = cfg.epsilon_proj_factor * h^2

    if cfg.use_projection_regularization is True, otherwise 0.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    cfg:
        LakeConfig instance.

    epsilon_proj:
        Optional projection regularization.

    Returns
    -------
    data:
        ProjectionData object.
    """
    if epsilon_proj is None:
        if cfg.use_projection_regularization:
            epsilon_proj = cfg.epsilon_proj_factor * grid.h**2
        else:
            epsilon_proj = 0.0

    projection_unknowns, idx_proj = build_projection_unknowns(grid)

    b_proj = build_projection_coefficient(
        grid=grid,
        epsilon_proj=epsilon_proj,
    )

    data = assemble_degenerate_projection_matrix(
        grid=grid,
        b_proj=b_proj,
        projection_unknowns=projection_unknowns,
        idx_proj=idx_proj,
        epsilon_proj=epsilon_proj,
    )

    return data


def print_projection_matrix_summary(
    data: ProjectionData,
) -> None:
    """
    Print a compact summary of the projection matrix.

    Parameters
    ----------
    data:
        ProjectionData object.
    """
    print("=" * 70)
    print("Degenerate projection matrix summary")
    print("=" * 70)
    print("Operator: A_proj = -div(b_proj grad)")
    print("Face averaging: arithmetic average of b_proj")
    print("Boundary treatment: first-order masked-grid Dirichlet flux")
    print(f"epsilon_proj      = {data.epsilon_proj:.8e}")
    print(f"shape             = {data.A_proj.shape}")
    print(f"number unknowns   = {data.num_unknowns}")
    print(f"nnz               = {data.nnz}")
    print(f"symmetry defect   = {data.symmetry_defect:.8e}")

    if data.epsilon_proj > 0.0:
        print("Note: b_proj = b + epsilon_proj.")
        print("      Projection is stabilized but approximate.")
    else:
        print("Note: b_proj = b.")
        print("      Degenerate coefficient used without regularization.")
        print("      Conditioning may be poor.")

    print("=" * 70)


# =============================================================================
# Weighted projection
# =============================================================================

def weighted_projection(
    u_star: np.ndarray,
    grid: LakeGrid,
    projection_data: ProjectionData,
    cfg: LakeConfig,
) -> ProjectionResult:
    """
    Apply the weighted projection to a provisional velocity u_star.

    Ideal continuous equation:

        div(b grad phi) = div(b u_star),

        u_projected = u_star - grad phi.

    SPD discrete solve:

        A_proj phi = -div(b u_star),

    where:

        A_proj = -div(b_proj grad).

    If b_proj != b, the projection is stabilized but approximate.

    Important
    ---------
    This is not an exact discrete weighted Helmholtz projection.

    The sparse operator A_proj is assembled in divergence form, but the
    correction

        u_projected = u_star - grad(phi)

    uses the prototype finite-difference gradient grad_scalar from
    operators.py.

    Therefore, exact discrete cancellation of div(b u_projected) is not
    guaranteed, even when epsilon_proj = 0.

    The quantities projection_defect_after and
    relative_projection_defect_after must be monitored.

    Parameters
    ----------
    u_star:
        Provisional velocity.

    grid:
        LakeGrid instance.

    projection_data:
        ProjectionData object.

    cfg:
        LakeConfig instance.

    Returns
    -------
    result:
        ProjectionResult object.
    """
    A = projection_data.A_proj

    div_before = div_bu(
        u_star,
        grid.b,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    rhs_grid = -div_before
    rhs_vec = rhs_grid[projection_data.projection_unknowns]

    phi_vec, cg_info, cg_iterations = cg_solve_compat(
        A=A,
        rhs=rhs_vec,
        rtol=cfg.cg_rtol,
        atol=cfg.cg_atol,
        maxiter=cfg.cg_maxiter,
    )

    phi = np.zeros_like(grid.b)
    phi[projection_data.projection_unknowns] = phi_vec
    phi[grid.boundary_band_topological] = 0.0
    phi[~grid.mask] = 0.0

    # Important:
    # grad_scalar is a lightweight prototype gradient based on numpy.gradient.
    # It is not the exact algebraic adjoint of the sparse divergence-form
    # operator used to assemble A_proj. Therefore this projection is only
    # approximately compatible at the discrete level.
    grad_phi = grad_scalar(
        phi,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    u_projected = u_star - grad_phi
    u_projected[~grid.mask] = 0.0

    div_after = div_bu(
        u_projected,
        grid.b,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    projection_defect_before = masked_l2_norm(
        div_before,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    projection_defect_after = masked_l2_norm(
        div_after,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    bu_projected_l2 = ordinary_l2_vector_norm(
        grid.b[..., None] * u_projected,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    relative_projection_defect_after = projection_defect_after / max(
        bu_projected_l2 / grid.h,
        1.0e-14,
    )

    residual_vec = A @ phi_vec - rhs_vec

    solver_residual = float(np.linalg.norm(residual_vec))
    rhs_norm = float(np.linalg.norm(rhs_vec))
    relative_solver_residual = solver_residual / max(rhs_norm, 1.0e-14)

    info = {
        "epsilon_proj": projection_data.epsilon_proj,
        "cg_info": cg_info,
        "cg_iterations": cg_iterations,
        "solver_residual": solver_residual,
        "relative_solver_residual": relative_solver_residual,
        "projection_defect_before": projection_defect_before,
        "projection_defect_after": projection_defect_after,
        "relative_projection_defect_after": relative_projection_defect_after,
        "rhs_norm": rhs_norm,
    }

    return ProjectionResult(
        u_projected=u_projected,
        phi=phi,
        info=info,
    )


# =============================================================================
# Projection diagnostics
# =============================================================================

def print_projection_diagnostics(
    result: ProjectionResult,
) -> None:
    """
    Print projection diagnostics.

    Parameters
    ----------
    result:
        ProjectionResult object.
    """
    info = result.info

    print("=" * 70)
    print("Weighted projection diagnostics")
    print("=" * 70)

    print(f"epsilon_proj                         = {info['epsilon_proj']:.8e}")
    print(f"CG info                              = {info['cg_info']}")
    print(f"CG iterations                        = {info['cg_iterations']}")
    print(f"solver residual                      = {info['solver_residual']:.8e}")
    print(f"relative solver residual             = {info['relative_solver_residual']:.8e}")
    print(f"projection defect before             = {info['projection_defect_before']:.8e}")
    print(f"projection defect after              = {info['projection_defect_after']:.8e}")
    print(f"relative projection defect after     = {info['relative_projection_defect_after']:.8e}")

    if info["epsilon_proj"] > 0.0:
        print()
        print("Warning:")
        print("  This projection used b_proj = b + epsilon_proj.")
        print("  It is stabilized but not an exact weighted projection.")
        print("  In addition, the correction gradient is not algebraically paired")
        print("  with A_proj, so div(bu) cancellation is only approximate.")
        print("  The projection defect must be tracked.")
    else:
        print()
        print("Note:")
        print("  This projection used b_proj = b.")
        print("  However, because the correction gradient is prototype-based,")
        print("  exact discrete div(bu) cancellation is still not guaranteed.")

    print("=" * 70)

""""

Ideal continuous projection:

    div(b grad phi) = div(b u_star),

    u_projected = u_star - grad phi.

If solved exactly with compatible continuous operators, this gives:

    div(b u_projected) = 0.

However, because b -> 0 near the boundary, the operator div(b grad phi) is
degenerate and can be ill-conditioned.

For numerical stability, this prototype allows a stabilized coefficient:

    b_proj = b + epsilon_proj.

Then the projection solve becomes:

    div(b_proj grad phi) = div(b u_star).

This is no longer an exact weighted Helmholtz projection, so the resulting
defect div(b u_projected) must be monitored.

Numerical SPD form:

    A_proj phi = rhs,

where:

    A_proj = -div(b_proj grad),
    rhs    = -div(b u_star).
"""