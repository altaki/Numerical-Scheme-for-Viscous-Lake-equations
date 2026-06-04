"""
Finite-difference operators and initial-field construction.

This module- perpendicular gradient,This module contains the basic discrete operators used throughout the project:
- scalar curl,
- vector divergence,
- weighted divergence div(bu),
- masked L2 and Linf norms,
- weighted L2_b velocity norm,
- weighted kinetic energy,
- construction of u0 and omega0.

Mathematical conventions used everywhere:

    grad_perp(psi) = (-partial_y psi, partial_x psi)

    curl(u) = partial_x u_y - partial_y u_x

    u = (1/b) grad_perp(psi)

    omega = curl(u) / b

    div((1/b) grad psi) = b omega

Important numerical note
------------------------
The differential operators in this module are lightweight Cartesian
finite-difference prototypes based on numpy.gradient followed by post-masking.

They are useful for initialization, diagnostics, plotting, and prototype
consistency checks. They are not boundary-fitted masked stencils and should not
be interpreted as research-grade operators near the degenerate boundary.

The boundary-sensitive sparse operators are assembled separately in:

    elliptic.py    for the singular inviscid operator -div((1/b)grad),
    projection.py  for the degenerate projection operator -div(b grad).
"""

from dataclasses import dataclass
from typing import Dict, Any

import numpy as np

from .grid import LakeGrid, support_to_boundary_distance


@dataclass
class LakeInitialFields:
    """
    Container for initial fields.

    Attributes
    ----------
    psi0:
        Initial stream function.

    u0:
        Initial velocity field, defined by:

            u0 = (1/b) grad_perp psi0.

    omega0:
        Local finite-difference prototype vorticity, defined by:

            omega0 = div((1/b) grad psi0) / b.

        This is a local prototype. The sparse elliptic reconstruction is
        assembled and tested independently in elliptic.py.

    diagnostics:
        Dictionary containing initial residual diagnostics.
    """

    psi0: np.ndarray
    u0: np.ndarray
    omega0: np.ndarray
    diagnostics: Dict[str, Any]


# =============================================================================
# Norms
# =============================================================================

def masked_l2_norm(
    f: np.ndarray,
    mask: np.ndarray,
    dx: float,
    dy: float,
) -> float:
    """
    Compute the masked discrete L2 norm of a scalar field.

    Definition:

        ||f||_L2 = sqrt(sum_mask |f|^2 dx dy).

    Parameters
    ----------
    f:
        Scalar field.

    mask:
        Active grid mask.

    dx, dy:
        Grid spacings.

    Returns
    -------
    norm:
        Masked L2 norm.
    """
    return float(np.sqrt(np.sum(np.abs(f[mask]) ** 2) * dx * dy))


def masked_linf_norm(
    f: np.ndarray,
    mask: np.ndarray,
) -> float:
    """
    Compute the masked Linf norm of a scalar field.

    Definition:

        ||f||_Linf = max_mask |f|.

    Parameters
    ----------
    f:
        Scalar field.

    mask:
        Active grid mask.

    Returns
    -------
    norm:
        Masked Linf norm.
    """
    if not np.any(mask):
        return 0.0

    return float(np.max(np.abs(f[mask])))


def ordinary_l2_vector_norm(
    v: np.ndarray,
    mask: np.ndarray,
    dx: float,
    dy: float,
) -> float:
    """
    Compute the ordinary masked L2 norm of a vector field.

    Definition:

        ||v||_L2 = sqrt(sum_mask |v|^2 dx dy).

    Parameters
    ----------
    v:
        Vector field with shape (..., 2) or more generally (..., k).

    mask:
        Active grid mask.

    dx, dy:
        Grid spacings.

    Returns
    -------
    norm:
        Ordinary masked L2 norm of the vector field.
    """
    mag2 = np.sum(v**2, axis=-1)
    return float(np.sqrt(np.sum(mag2[mask]) * dx * dy))


def weighted_l2b_velocity_norm(
    u: np.ndarray,
    b: np.ndarray,
    mask: np.ndarray,
    dx: float,
    dy: float,
) -> float:
    """
    Compute the weighted velocity norm ||u||_{L^2_b}.

    Definition:

        ||u||_{L^2_b}
        =
        sqrt(sum_mask b |u|^2 dx dy).

    Parameters
    ----------
    u:
        Velocity field with shape (..., 2).

    b:
        Bathymetry.

    mask:
        Active grid mask.

    dx, dy:
        Grid spacings.

    Returns
    -------
    norm:
        Weighted L2_b velocity norm.
    """
    mag2 = u[..., 0] ** 2 + u[..., 1] ** 2

    return float(
        np.sqrt(
            np.sum((b * mag2)[mask]) * dx * dy
        )
    )


def weighted_kinetic_energy(
    u: np.ndarray,
    b: np.ndarray,
    mask: np.ndarray,
    dx: float,
    dy: float,
) -> float:
    """
    Compute the weighted kinetic energy.

    Definition:

        E = 1/2 int b |u|^2 dx.

    Discrete version:

        E_h = 1/2 sum_mask b |u|^2 dx dy.

    Parameters
    ----------
    u:
        Velocity field.

    b:
        Bathymetry.

    mask:
        Active mask.

    dx, dy:
        Grid spacings.

    Returns
    -------
    energy:
        Weighted kinetic energy.
    """
    mag2 = u[..., 0] ** 2 + u[..., 1] ** 2

    return float(
        0.5 * np.sum((b * mag2)[mask]) * dx * dy
    )


# =============================================================================
# Differential operators
# =============================================================================

def grad_scalar(
    f: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Compute a prototype finite-difference gradient of a scalar field.

    The derivative is computed on the full Cartesian grid using numpy.gradient,
    then values outside the active mask are set to zero.

    Important
    ---------
    This is a full-grid finite-difference prototype with post-masking.
    It is not a boundary-fitted masked stencil. Near the discrete boundary,
    centered differences may involve values outside the active mask.

    This is acceptable for the current initial-field prototype because psi0 is
    compactly supported away from the boundary. It should not be interpreted as
    a research-grade boundary operator.

    Parameters
    ----------
    f:
        Scalar field.

    dx, dy:
        Grid spacings.

    mask:
        Active grid mask.

    Returns
    -------
    grad_f:
        Vector field with components:

            grad_f[..., 0] = partial_x f,
            grad_f[..., 1] = partial_y f.
    """
    df_dy, df_dx = np.gradient(f, dy, dx)

    df_dx = np.asarray(df_dx)
    df_dy = np.asarray(df_dy)

    df_dx[~mask] = 0.0
    df_dy[~mask] = 0.0

    return np.stack([df_dx, df_dy], axis=-1)


def grad_perp_scalar(
    f: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Compute the perpendicular gradient using the fixed convention.

    Convention:

        grad_perp f = (-partial_y f, partial_x f).

    Parameters
    ----------
    f:
        Scalar field.

    dx, dy:
        Grid spacings.

    mask:
        Active mask.

    Returns
    -------
    grad_perp_f:
        Vector field.
    """
    grad_f = grad_scalar(f, dx, dy, mask)

    partial_x = grad_f[..., 0]
    partial_y = grad_f[..., 1]

    return np.stack([-partial_y, partial_x], axis=-1)


def curl_vector(
    u: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Compute the scalar curl of a vector field.

    Convention:

        curl(u) = partial_x u_y - partial_y u_x.

    Important
    ---------
    This is a full-grid finite-difference prototype with post-masking.
    It is not a boundary-fitted masked stencil. Near the discrete boundary,
    finite-difference stencils may involve values outside the active mask.

    Parameters
    ----------
    u:
        Vector field with components u[..., 0], u[..., 1].

    dx, dy:
        Grid spacings.

    mask:
        Active mask.

    Returns
    -------
    curl_u:
        Scalar curl field.
    """
    ux = u[..., 0]
    uy = u[..., 1]

    duy_dx = np.gradient(uy, dx, axis=1)
    dux_dy = np.gradient(ux, dy, axis=0)

    curl_u = duy_dx - dux_dy
    curl_u[~mask] = 0.0

    return curl_u


def div_vector(
    q: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Compute the divergence of a vector field.

    Definition:

        div(q) = partial_x q_x + partial_y q_y.

    Important
    ---------
    This is a full-grid finite-difference prototype with post-masking.
    It is not a boundary-fitted masked stencil. Near the discrete boundary,
    finite-difference stencils may involve values outside the active mask.

    Parameters
    ----------
    q:
        Vector field.

    dx, dy:
        Grid spacings.

    mask:
        Active mask.

    Returns
    -------
    div_q:
        Scalar divergence field.
    """
    qx = q[..., 0]
    qy = q[..., 1]

    dqx_dx = np.gradient(qx, dx, axis=1)
    dqy_dy = np.gradient(qy, dy, axis=0)

    div_q = dqx_dx + dqy_dy
    div_q[~mask] = 0.0

    return div_q


def div_bu(
    u: np.ndarray,
    b: np.ndarray,
    dx: float,
    dy: float,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Compute the weighted divergence div(bu).

    Definition:

        div(bu) = partial_x(b u_x) + partial_y(b u_y).

    This is the discrete prototype of the weighted incompressibility constraint:

        div(bu) = 0.

    Parameters
    ----------
    u:
        Velocity field.

    b:
        Bathymetry.

    dx, dy:
        Grid spacings.

    mask:
        Active mask.

    Returns
    -------
    div_bu_field:
        Scalar field div(bu).
    """
    return div_vector(b[..., None] * u, dx, dy, mask)


# =============================================================================
# Plot helper
# =============================================================================

def mask_to_nan(
    arr: np.ndarray,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Return an array equal to arr inside the mask and NaN outside.

    This helper works for:

    - scalar fields with shape (Ny, Nx);
    - vector fields with shape (Ny, Nx, k).

    Parameters
    ----------
    arr:
        Scalar or vector array.

    mask:
        Boolean mask with shape (Ny, Nx).

    Returns
    -------
    masked_arr:
        Array equal to arr inside the mask and NaN outside.

    Raises
    ------
    ValueError:
        If arr and mask have incompatible dimensions.
    """
    if arr.ndim == mask.ndim:
        return np.where(mask, arr, np.nan)

    if arr.ndim == mask.ndim + 1:
        return np.where(mask[..., None], arr, np.nan)

    raise ValueError(
        "mask_to_nan: incompatible array and mask dimensions. "
        f"arr.ndim={arr.ndim}, mask.ndim={mask.ndim}."
    )


# =============================================================================
# Initial fields
# =============================================================================

def build_initial_velocity(
    grid: LakeGrid,
) -> np.ndarray:
    """
    Build the initial velocity:

        u0 = (1/b) grad_perp psi0.

    In the code, b_safe is used defensively for division. This is harmless here
    because psi0 is compactly supported away from the boundary, so grad(psi0)
    vanishes near the region where b is close to zero.

    This should not be interpreted as a general-purpose regularization of the
    singular inviscid reconstruction operator.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    Returns
    -------
    u0:
        Initial velocity field.
    """
    grad_perp_psi0 = grad_perp_scalar(
        grid.psi0,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    u0 = grad_perp_psi0 / grid.b_safe[..., None]
    u0[~grid.mask] = 0.0

    return u0


def build_initial_vorticity_prototype(
    grid: LakeGrid,
) -> np.ndarray:
    """
    Build the local finite-difference prototype vorticity omega0.

    Definition:

        omega0 = div((1/b) grad psi0) / b.

    This is a local finite-difference prototype using numpy.gradient-based
    operators. It is useful for initialization and consistency diagnostics.

    It is not the authoritative sparse elliptic discretization. The singular
    elliptic reconstruction is assembled separately in elliptic.py through the
    divergence-form operator:

        A_sing = -div((1/b) grad).

    Parameters
    ----------
    grid:
        LakeGrid instance.

    Returns
    -------
    omega0:
        Prototype vorticity.
    """
    grad_psi0 = grad_scalar(
        grid.psi0,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    q0 = grad_psi0 / grid.b_safe[..., None]

    omega0 = div_vector(
        q0,
        grid.dx,
        grid.dy,
        grid.mask,
    ) / grid.b_safe

    omega0[~grid.mask] = 0.0

    return omega0


def build_initial_fields(
    grid: LakeGrid,
) -> LakeInitialFields:
    """
    Construct psi0, u0, omega0 and diagnostic residuals.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    Returns
    -------
    fields:
        LakeInitialFields object.
    """
    psi0 = grid.psi0.copy()

    u0 = build_initial_velocity(grid)

    omega0 = build_initial_vorticity_prototype(grid)

    div_bu0 = div_bu(
        u0,
        grid.b,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    curl_u0 = curl_vector(
        u0,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    curl_over_b_u0 = curl_u0 / grid.b_safe

    res_vorticity0 = curl_over_b_u0 - omega0

    div_bu0_l2 = masked_l2_norm(
        div_bu0,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    vorticity_residual_l2 = masked_l2_norm(
        res_vorticity0,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    omega0_l2 = masked_l2_norm(
        omega0,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    relative_vorticity_residual = vorticity_residual_l2 / max(
        omega0_l2,
        1.0e-14,
    )

    bu0_l2 = ordinary_l2_vector_norm(
        grid.b[..., None] * u0,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    relative_div_bu0 = div_bu0_l2 / max(
        bu0_l2 / grid.h,
        1.0e-14,
    )

    diagnostics = {
        "div_bu0_l2": div_bu0_l2,
        "relative_div_bu0": relative_div_bu0,
        "vorticity_residual_l2": vorticity_residual_l2,
        "relative_vorticity_residual": relative_vorticity_residual,
        "omega0_l2": omega0_l2,
        "support_distance_psi0": support_to_boundary_distance(
            psi0,
            grid,
            threshold=1.0e-6,
        ),
        "support_distance_omega0": support_to_boundary_distance(
            omega0,
            grid,
            threshold=1.0e-3,
        ),
        "weighted_energy_u0": weighted_kinetic_energy(
            u0,
            grid.b,
            grid.mask,
            grid.dx,
            grid.dy,
        ),
        "weighted_l2b_norm_u0": weighted_l2b_velocity_norm(
            u0,
            grid.b,
            grid.mask,
            grid.dx,
            grid.dy,
        ),
    }

    return LakeInitialFields(
        psi0=psi0,
        u0=u0,
        omega0=omega0,
        diagnostics=diagnostics,
    )


def print_initial_diagnostics(
    fields: LakeInitialFields,
) -> None:
    """
    Print initial-field diagnostics.

    Parameters
    ----------
    fields:
        LakeInitialFields instance.
    """
    d = fields.diagnostics

    print("=" * 70)
    print("Initial-field diagnostics")
    print("=" * 70)

    print(f"||div(bu0)||_L2                       = {d['div_bu0_l2']:.8e}")
    print(f"relative div(bu0)                     = {d['relative_div_bu0']:.8e}")
    print(f"||curl(u0)/b - omega0||_L2            = {d['vorticity_residual_l2']:.8e}")
    print(f"relative vorticity residual           = {d['relative_vorticity_residual']:.8e}")
    print(f"||omega0||_L2                         = {d['omega0_l2']:.8e}")
    print(f"support distance psi0                 = {d['support_distance_psi0']:.8e}")
    print(f"support distance omega0               = {d['support_distance_omega0']:.8e}")
    print(f"weighted kinetic energy of u0          = {d['weighted_energy_u0']:.8e}")
    print(f"||u0||_L2_b                            = {d['weighted_l2b_norm_u0']:.8e}")

    print("=" * 70)


