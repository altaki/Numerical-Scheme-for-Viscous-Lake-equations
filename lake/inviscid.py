"""
Inviscid lake-equation solver in vorticity formulation.

This module implements the inviscid lake equations:

    partial_t omega + u · grad omega = 0,

with velocity reconstructed from the singular elliptic problem:

    div((1/b) grad psi) = b omega,

and

    u = (1/b) grad_perp psi.

The reconstruction is performed by elliptic.reconstruct_from_vorticity.

The time transport is done by a semi-Lagrangian prototype scheme.

Important:
    This solver is for the inviscid problem only.
    It must not be used for the viscous equations.

Mathematical conventions:

    grad_perp(psi) = (-partial_y psi, partial_x psi)

    curl(u) = partial_x u_y - partial_y u_x

    omega = curl(u)/b

    div((1/b) grad psi) = b omega
"""

from dataclasses import dataclass
from typing import Dict, Any

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .config import LakeConfig
from .grid import LakeGrid, support_to_boundary_distance
from .elliptic import SingularEllipticData, reconstruct_from_vorticity
from .operators import (
    curl_vector,
    div_bu,
    masked_l2_norm,
    masked_linf_norm,
    ordinary_l2_vector_norm,
)


@dataclass
class InviscidResult:
    """
    Result of the inviscid vorticity solver.

    Attributes
    ----------
    omega_final:
        Final vorticity field.

    psi_final:
        Final stream function.

    u_final:
        Final reconstructed velocity field.

    diagnostics:
        Dictionary of time histories.
    """

    omega_final: np.ndarray
    psi_final: np.ndarray
    u_final: np.ndarray
    diagnostics: Dict[str, np.ndarray]


# =============================================================================
# Semi-Lagrangian transport
# =============================================================================

def semi_lagrangian_step(
    omega: np.ndarray,
    u: np.ndarray,
    grid: LakeGrid,
    dt: float,
) -> np.ndarray:
    """
    Perform one semi-Lagrangian step for:

        partial_t omega + u · grad omega = 0.

    Backward characteristics:

        X_back = X - dt u_x,
        Y_back = Y - dt u_y.

    Interpolation is done with scipy.interpolate.RegularGridInterpolator.

    Important array convention:
        omega has shape (len(y), len(x)),
        so the interpolator grid order is (y, x).

    Parameters
    ----------
    omega:
        Current vorticity field.

    u:
        Current velocity field.

    grid:
        LakeGrid instance.

    dt:
        Time step.

    Returns
    -------
    omega_new:
        Transported vorticity.
    """
    X_back = grid.X - dt * u[..., 0]
    Y_back = grid.Y - dt * u[..., 1]

    interpolator = RegularGridInterpolator(
        (grid.y, grid.x),
        omega,
        bounds_error=False,
        fill_value=0.0,
    )

    points = np.stack(
        [Y_back.ravel(), X_back.ravel()],
        axis=-1,
    )

    omega_new = interpolator(points).reshape(omega.shape)

    omega_new[~grid.mask] = 0.0

    return omega_new


# =============================================================================
# Diagnostics
# =============================================================================

def compute_inviscid_diagnostics(
    omega: np.ndarray,
    psi: np.ndarray,
    u: np.ndarray,
    reconstruction_info: Dict[str, Any],
    grid: LakeGrid,
    time: float,
    dt: float,
) -> Dict[str, float]:
    """
    Compute one snapshot of inviscid diagnostics.

    Parameters
    ----------
    omega:
        Vorticity field.

    psi:
        Stream function.

    u:
        Velocity field.

    reconstruction_info:
        Information dictionary returned by reconstruct_from_vorticity.

    grid:
        LakeGrid instance.

    time:
        Current time.

    dt:
        Time step used to reach the next state.

    Returns
    -------
    diag:
        Dictionary of scalar diagnostics.
    """
    speed = np.sqrt(u[..., 0] ** 2 + u[..., 1] ** 2)
    max_speed = float(np.max(speed[grid.mask])) if np.any(grid.mask) else 0.0

    div_val = div_bu(
        u,
        grid.b,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    div_l2 = masked_l2_norm(
        div_val,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    bu_l2 = ordinary_l2_vector_norm(
        grid.b[..., None] * u,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    relative_div = div_l2 / max(bu_l2 / grid.h, 1.0e-14)

    curl_u = curl_vector(
        u,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    curl_over_b = curl_u / grid.b_safe

    vorticity_residual = curl_over_b - omega

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

    support_distance = support_to_boundary_distance(
        omega,
        grid,
        threshold=1.0e-3,
    )

    diag = {
        "time": float(time),
        "dt": float(dt),
        "min_omega": float(np.min(omega[grid.mask])) if np.any(grid.mask) else np.nan,
        "max_omega": float(np.max(omega[grid.mask])) if np.any(grid.mask) else np.nan,
        "omega_linf": masked_linf_norm(omega, grid.mask),
        "omega_l2": omega_l2,
        "support_distance": support_distance,
        "divbu_l2": div_l2,
        "relative_divbu": relative_div,
        "vorticity_residual_l2": vorticity_residual_l2,
        "relative_vorticity_residual": relative_vorticity_residual,
        "elliptic_residual": float(reconstruction_info["elliptic_residual"]),
        "relative_elliptic_residual": float(
            reconstruction_info["relative_elliptic_residual"]
        ),
        "cg_iterations": int(reconstruction_info["cg_iterations"]),
        "max_speed": max_speed,
        "cfl_ratio": max_speed * dt / grid.h if grid.h > 0 else np.nan,
    }

    return diag


def _append_diag(history: Dict[str, list], diag: Dict[str, float]) -> None:
    """
    Append scalar diagnostic values to a history dictionary.

    Parameters
    ----------
    history:
        Dictionary of lists.

    diag:
        Dictionary of scalar diagnostics.
    """
    for key, value in diag.items():
        history.setdefault(key, []).append(value)


def _finalize_diag(history: Dict[str, list]) -> Dict[str, np.ndarray]:
    """
    Convert diagnostic history lists to numpy arrays.

    Parameters
    ----------
    history:
        Dictionary of lists.

    Returns
    -------
    diagnostics:
        Dictionary of numpy arrays.
    """
    return {key: np.asarray(value) for key, value in history.items()}


# =============================================================================
# Main inviscid solver
# =============================================================================

def run_inviscid_solver(
    omega_initial: np.ndarray,
    grid: LakeGrid,
    elliptic_data: SingularEllipticData,
    cfg: LakeConfig,
    T: float | None = None,
    CFL: float | None = None,
    dt_max: float | None = None,
) -> InviscidResult:
    """
    Run the inviscid lake solver in vorticity formulation.

    Equation:

        partial_t omega + u · grad omega = 0.

    Velocity reconstruction:

        div((1/b) grad psi) = b omega,
        u = (1/b) grad_perp psi.

    Parameters
    ----------
    omega_initial:
        Initial vorticity.

    grid:
        LakeGrid instance.

    elliptic_data:
        Singular elliptic matrix data.

    cfg:
        LakeConfig instance.

    T:
        Final time. If None, uses cfg.T_inviscid.

    CFL:
        CFL-like factor. If None, uses cfg.CFL_inviscid.

    dt_max:
        Maximum timestep. If None, uses cfg.dt_max_inviscid.

    Returns
    -------
    result:
        InviscidResult object.
    """
    if T is None:
        T = cfg.T_inviscid

    if CFL is None:
        CFL = cfg.CFL_inviscid

    if dt_max is None:
        dt_max = cfg.dt_max_inviscid

    omega = omega_initial.copy()
    omega[~grid.mask] = 0.0

    diagnostics_history: Dict[str, list] = {}

    t = 0.0

    while t < T - 1.0e-15:
        reconstruction = reconstruct_from_vorticity(
            omega=omega,
            grid=grid,
            elliptic_data=elliptic_data,
            cfg=cfg,
        )

        psi = reconstruction.psi
        u = reconstruction.u
        info = reconstruction.info

        speed = np.sqrt(u[..., 0] ** 2 + u[..., 1] ** 2)
        max_speed = float(np.max(speed[grid.mask])) if np.any(grid.mask) else 0.0

        dt = min(
            dt_max,
            CFL * grid.h / max(max_speed, 1.0e-14),
            T - t,
        )

        diag = compute_inviscid_diagnostics(
            omega=omega,
            psi=psi,
            u=u,
            reconstruction_info=info,
            grid=grid,
            time=t,
            dt=dt,
        )

        _append_diag(diagnostics_history, diag)

        omega = semi_lagrangian_step(
            omega=omega,
            u=u,
            grid=grid,
            dt=dt,
        )

        omega[~grid.mask] = 0.0

        t += dt

    # Final reconstruction at time T
    reconstruction_final = reconstruct_from_vorticity(
        omega=omega,
        grid=grid,
        elliptic_data=elliptic_data,
        cfg=cfg,
    )

    final_diag = compute_inviscid_diagnostics(
        omega=omega,
        psi=reconstruction_final.psi,
        u=reconstruction_final.u,
        reconstruction_info=reconstruction_final.info,
        grid=grid,
        time=t,
        dt=0.0,
    )

    _append_diag(diagnostics_history, final_diag)

    diagnostics = _finalize_diag(diagnostics_history)

    return InviscidResult(
        omega_final=omega,
        psi_final=reconstruction_final.psi,
        u_final=reconstruction_final.u,
        diagnostics=diagnostics,
    )


# =============================================================================
# Printing and quick diagnostics
# =============================================================================

def print_inviscid_summary(
    result: InviscidResult,
) -> None:
    """
    Print a compact summary of an inviscid simulation.

    Parameters
    ----------
    result:
        InviscidResult object.
    """
    d = result.diagnostics

    print("=" * 70)
    print("Inviscid solver summary")
    print("=" * 70)

    if "time" in d and len(d["time"]) > 0:
        print(f"final time                         = {d['time'][-1]:.8e}")
        print(f"number of stored steps             = {len(d['time'])}")

    if "min_omega" in d:
        print(f"final min omega                    = {d['min_omega'][-1]:.8e}")

    if "max_omega" in d:
        print(f"final max omega                    = {d['max_omega'][-1]:.8e}")

    if "omega_linf" in d:
        print(f"final ||omega||_Linf               = {d['omega_linf'][-1]:.8e}")

    if "support_distance" in d:
        print(f"final support distance             = {d['support_distance'][-1]:.8e}")

    if "divbu_l2" in d:
        print(f"final ||div(bu)||_L2               = {d['divbu_l2'][-1]:.8e}")

    if "relative_divbu" in d:
        print(f"final relative div(bu)             = {d['relative_divbu'][-1]:.8e}")

    if "vorticity_residual_l2" in d:
        print(f"final ||curl(u)/b - omega||_L2     = {d['vorticity_residual_l2'][-1]:.8e}")

    if "elliptic_residual" in d:
        print(f"final elliptic residual            = {d['elliptic_residual'][-1]:.8e}")

    if "cg_iterations" in d:
        print(f"average CG iterations              = {np.mean(d['cg_iterations']):.3f}")

    if "cfl_ratio" in d:
        print(f"max CFL ratio                      = {np.max(d['cfl_ratio']):.8e}")

    print("=" * 70)