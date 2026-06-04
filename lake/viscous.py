"""
Prototype viscous lake-equation solver in velocity formulation.

This module implements a masked Cartesian predictor-projection prototype for
the viscous lake equations.

Continuous model:

    partial_t(b u_mu)
    + div(b u_mu tensor u_mu)
    - 2 mu div(b D(u_mu) + b div(u_mu) I)
    + b grad(p_mu)
    = 0,

with constraint:

    div(b u_mu) = 0.

The continuous boundary conditions are Navier-type:

    b u_mu · n = 0,

and

    2 b (D(u_mu) · n + div(u_mu)n) · tau
    + eta_mu b(u_mu · tau) = 0.

Important:
    This prototype does NOT exactly impose the continuous Navier boundary
    condition. It uses a masked Cartesian approximation and monitors the
    resulting weighted-divergence/projection defects.

Numerical method:

    1. Explicit/semi-explicit predictor:

        u_star = u
                 - dt * div(b u tensor u) / b_safe
                 + dt * mu * viscous_operator(u) / b_safe

       where viscous_operator approximates:

        2 div(bD(u) + b div(u) I).

    2. Degenerate weighted projection:

        div(b grad phi) = div(b u_star),

       implemented through projection.weighted_projection.

    3. Velocity update:

        u_new = u_star - grad(phi).

The solver uses velocity formulation only.
It does not use the viscous vorticity equation.
"""

from dataclasses import dataclass
from typing import Dict, Any, Optional

import numpy as np

from .config import LakeConfig
from .grid import LakeGrid
from .operators import (
    div_bu,
    masked_l2_norm,
    weighted_kinetic_energy,
    ordinary_l2_vector_norm,
)
from .projection import ProjectionData, weighted_projection


@dataclass
class ViscousResult:
    """
    Result of the viscous velocity solver.

    Attributes
    ----------
    u_final:
        Final viscous velocity.

    p_final:
        Final pressure/projection-potential approximation.

        In the projection method, the projection potential phi satisfies
        u_new = u_star - grad(phi). A pressure-like quantity is approximated
        by phi/dt.

    diagnostics:
        Dictionary of time histories.
    """

    u_final: np.ndarray
    p_final: np.ndarray
    diagnostics: Dict[str, np.ndarray]


# =============================================================================
# Basic differential terms for viscous velocity solver
# =============================================================================

def convective_term(
    u: np.ndarray,
    b: np.ndarray,
    grid: LakeGrid,
) -> np.ndarray:
    """
    Compute a prototype conservative approximation of div(b u tensor u).

    For u=(u_x,u_y), define:

        F_xx = b u_x u_x
        F_xy = b u_x u_y
        F_yx = b u_y u_x
        F_yy = b u_y u_y

    Then:

        conv_x = partial_x F_xx + partial_y F_xy
        conv_y = partial_x F_yx + partial_y F_yy.

    Parameters
    ----------
    u:
        Velocity field.

    b:
        Bathymetry.

    grid:
        LakeGrid instance.

    Returns
    -------
    conv:
        Vector field approximating div(b u tensor u).
    """
    ux = u[..., 0]
    uy = u[..., 1]

    Fxx = b * ux * ux
    Fxy = b * ux * uy
    Fyx = b * uy * ux
    Fyy = b * uy * uy

    dFxx_dy, dFxx_dx = np.gradient(Fxx, grid.dy, grid.dx)
    dFxy_dy, dFxy_dx = np.gradient(Fxy, grid.dy, grid.dx)
    dFyx_dy, dFyx_dx = np.gradient(Fyx, grid.dy, grid.dx)
    dFyy_dy, dFyy_dx = np.gradient(Fyy, grid.dy, grid.dx)

    conv_x = dFxx_dx + dFxy_dy
    conv_y = dFyx_dx + dFyy_dy

    conv = np.stack([conv_x, conv_y], axis=-1)
    conv[~grid.mask] = 0.0

    return conv


def div_u(
    u: np.ndarray,
    grid: LakeGrid,
) -> np.ndarray:
    """
    Compute ordinary divergence div(u).

    Parameters
    ----------
    u:
        Velocity field.

    grid:
        LakeGrid instance.

    Returns
    -------
    divu:
        Scalar divergence field.
    """
    ux = u[..., 0]
    uy = u[..., 1]

    dux_dy, dux_dx = np.gradient(ux, grid.dy, grid.dx)
    duy_dy, duy_dx = np.gradient(uy, grid.dy, grid.dx)

    divu = dux_dx + duy_dy
    divu[~grid.mask] = 0.0

    return divu


def viscous_operator(
    u: np.ndarray,
    b: np.ndarray,
    grid: LakeGrid,
) -> np.ndarray:
    """
    Prototype approximation of the viscous operator:

        2 div(bD(u) + b div(u) I).

    With:

        D(u) = 0.5(grad u + grad u^T).

    Components:

        D11 = partial_x u_x
        D22 = partial_y u_y
        D12 = 0.5(partial_y u_x + partial_x u_y)

        divu = partial_x u_x + partial_y u_y

        T11 = b(D11 + divu)
        T12 = bD12
        T21 = bD12
        T22 = b(D22 + divu)

        visc_x = 2(partial_x T11 + partial_y T12)
        visc_y = 2(partial_x T21 + partial_y T22)

    Parameters
    ----------
    u:
        Velocity field.

    b:
        Bathymetry.

    grid:
        LakeGrid instance.

    Returns
    -------
    visc:
        Vector field approximating the viscous operator.
    """
    ux = u[..., 0]
    uy = u[..., 1]

    dux_dy, dux_dx = np.gradient(ux, grid.dy, grid.dx)
    duy_dy, duy_dx = np.gradient(uy, grid.dy, grid.dx)

    D11 = dux_dx
    D22 = duy_dy
    D12 = 0.5 * (dux_dy + duy_dx)

    divu = dux_dx + duy_dy

    T11 = b * (D11 + divu)
    T12 = b * D12
    T21 = b * D12
    T22 = b * (D22 + divu)

    dT11_dy, dT11_dx = np.gradient(T11, grid.dy, grid.dx)
    dT12_dy, dT12_dx = np.gradient(T12, grid.dy, grid.dx)
    dT21_dy, dT21_dx = np.gradient(T21, grid.dy, grid.dx)
    dT22_dy, dT22_dx = np.gradient(T22, grid.dy, grid.dx)

    visc_x = 2.0 * (dT11_dx + dT12_dy)
    visc_y = 2.0 * (dT21_dx + dT22_dy)

    visc = np.stack([visc_x, visc_y], axis=-1)
    visc[~grid.mask] = 0.0

    return visc


# =============================================================================
# Boundary diagnostic
# =============================================================================

def no_penetration_defect(
    u: np.ndarray,
    grid: LakeGrid,
) -> Dict[str, float]:
    """
    Approximate the no-penetration defect on the diagnostic boundary band.

    Continuous condition:

        b u · n = 0.

    For the unit disk, approximate:

        n = (x/r, y/r)

    on the geometric boundary band grid.boundary_band.

    Parameters
    ----------
    u:
        Velocity field.

    grid:
        LakeGrid instance.

    Returns
    -------
    defect:
        Dictionary with max and L2-like boundary-band defects.
    """
    band = grid.boundary_band

    if not np.any(band):
        return {"max_abs": np.nan, "l2": np.nan}

    nx = np.zeros_like(grid.r)
    ny = np.zeros_like(grid.r)

    safe_r = np.maximum(grid.r, 1.0e-14)
    nx[band] = grid.X[band] / safe_r[band]
    ny[band] = grid.Y[band] / safe_r[band]

    bu_dot_n = grid.b * (u[..., 0] * nx + u[..., 1] * ny)

    max_abs = float(np.max(np.abs(bu_dot_n[band])))

    # Boundary-band L2-like diagnostic, not a true boundary integral.
    l2 = float(np.sqrt(np.sum(bu_dot_n[band] ** 2) * grid.dx * grid.dy))

    return {
        "max_abs": max_abs,
        "l2": l2,
    }


# =============================================================================
# Diagnostics
# =============================================================================

def relative_div_bu(
    u: np.ndarray,
    grid: LakeGrid,
) -> float:
    """
    Compute relative weighted-divergence defect.

    Definition:

        ||div(bu)||_L2 / max(||bu||_L2/h, 1e-14).

    Parameters
    ----------
    u:
        Velocity field.

    grid:
        LakeGrid instance.

    Returns
    -------
    relative_defect:
        Relative divergence defect.
    """
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

    return div_l2 / max(bu_l2 / grid.h, 1.0e-14)


def compute_viscous_snapshot_diagnostics(
    u: np.ndarray,
    projection_info: Dict[str, Any],
    mu: float,
    dt: float,
    time: float,
    grid: LakeGrid,
) -> Dict[str, float]:
    """
    Compute diagnostics for one viscous time step.

    Parameters
    ----------
    u:
        Current velocity.

    projection_info:
        Diagnostics from weighted_projection.

    mu:
        Viscosity.

    dt:
        Time step.

    time:
        Current time.

    grid:
        LakeGrid instance.

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

    boundary_defect = no_penetration_defect(u, grid)

    diag = {
        "time": float(time),
        "dt": float(dt),
        "mu": float(mu),
        "weighted_energy": weighted_kinetic_energy(
            u,
            grid.b,
            grid.mask,
            grid.dx,
            grid.dy,
        ),
        "divbu_l2": float(div_l2),
        "relative_divbu": float(relative_div_bu(u, grid)),
        "projection_defect_before": float(
            projection_info.get("projection_defect_before", np.nan)
        ),
        "projection_defect_after": float(
            projection_info.get("projection_defect_after", np.nan)
        ),
        "projection_solver_residual": float(
            projection_info.get("solver_residual", np.nan)
        ),
        "relative_projection_solver_residual": float(
            projection_info.get("relative_solver_residual", np.nan)
        ),
        "projection_cg_iterations": int(
            projection_info.get("cg_iterations", -1)
        ),
        "max_speed": max_speed,
        "cfl_ratio": max_speed * dt / grid.h if grid.h > 0 else np.nan,
        "boundary_no_penetration_max": boundary_defect["max_abs"],
        "boundary_no_penetration_l2": boundary_defect["l2"],
    }

    return diag


def _append_diag(history: Dict[str, list], diag: Dict[str, float]) -> None:
    """
    Append a scalar diagnostic dictionary to a dictionary of lists.
    """
    for key, value in diag.items():
        history.setdefault(key, []).append(value)


def _finalize_diag(history: Dict[str, list]) -> Dict[str, np.ndarray]:
    """
    Convert lists in diagnostic history to numpy arrays.
    """
    return {key: np.asarray(value) for key, value in history.items()}


# =============================================================================
# Viscous solver
# =============================================================================

def viscous_predictor_explicit(
    u: np.ndarray,
    mu: float,
    dt: float,
    grid: LakeGrid,
) -> np.ndarray:
    """
    Explicit/semi-explicit predictor for the viscous velocity equation.

    Prototype formula:

        u_star = u
                 - dt * div(b u tensor u) / b_safe
                 + dt * mu * viscous_operator(u) / b_safe.

    Parameters
    ----------
    u:
        Current velocity.

    mu:
        Viscosity.

    dt:
        Time step.

    grid:
        LakeGrid instance.

    Returns
    -------
    u_star:
        Provisional velocity.
    """
    conv = convective_term(
        u,
        grid.b,
        grid,
    )

    visc = viscous_operator(
        u,
        grid.b,
        grid,
    )

    u_star = (
        u
        - dt * conv / grid.b_safe[..., None]
        + dt * mu * visc / grid.b_safe[..., None]
    )

    u_star[~grid.mask] = 0.0

    return u_star


def choose_viscous_timestep(
    u: np.ndarray,
    mu: float,
    grid: LakeGrid,
    cfg: LakeConfig,
    T_remaining: float,
    CFL: Optional[float] = None,
    dt_max: Optional[float] = None,
) -> float:
    """
    Choose a conservative time step for the explicit viscous prototype.

    Restrictions:

        dt_adv  = CFL h / max|u|
        dt_visc = 0.25 h^2 / mu
        dt      = min(dt_max, dt_adv, dt_visc, T_remaining)

    Parameters
    ----------
    u:
        Current velocity.

    mu:
        Viscosity.

    grid:
        LakeGrid instance.

    cfg:
        LakeConfig instance.

    T_remaining:
        Remaining time.

    CFL:
        Optional CFL factor.

    dt_max:
        Optional maximum timestep.

    Returns
    -------
    dt:
        Chosen timestep.
    """
    if CFL is None:
        CFL = cfg.CFL_viscous

    if dt_max is None:
        dt_max = cfg.dt_max_viscous

    speed = np.sqrt(u[..., 0] ** 2 + u[..., 1] ** 2)
    max_speed = float(np.max(speed[grid.mask])) if np.any(grid.mask) else 0.0

    dt_adv = CFL * grid.h / max(max_speed, 1.0e-14)
    dt_visc = 0.25 * grid.h * grid.h / max(mu, 1.0e-14)

    dt = min(
        dt_max,
        dt_adv,
        dt_visc,
        T_remaining,
    )

    return float(dt)


def run_viscous_solver(
    u_initial: np.ndarray,
    mu: float,
    grid: LakeGrid,
    projection_data: ProjectionData,
    cfg: LakeConfig,
    T: Optional[float] = None,
    CFL: Optional[float] = None,
    dt_max: Optional[float] = None,
    project_initial: bool = True,
) -> ViscousResult:
    """
    Run the viscous lake solver in velocity formulation.

    This is a prototype predictor-projection solver.

    Parameters
    ----------
    u_initial:
        Initial velocity.

    mu:
        Viscosity.

    grid:
        LakeGrid instance.

    projection_data:
        ProjectionData for the degenerate weighted projection.

    cfg:
        LakeConfig instance.

    T:
        Final time. If None, uses cfg.T_viscous.

    CFL:
        CFL factor. If None, uses cfg.CFL_viscous.

    dt_max:
        Maximum timestep. If None, uses cfg.dt_max_viscous.

    project_initial:
        If True, project the initial velocity before time stepping.

    Returns
    -------
    result:
        ViscousResult object.
    """
    if T is None:
        T = cfg.T_viscous

    if CFL is None:
        CFL = cfg.CFL_viscous

    if dt_max is None:
        dt_max = cfg.dt_max_viscous

    u = u_initial.copy()
    u[~grid.mask] = 0.0

    p_current = np.zeros_like(grid.b)

    diagnostics_history: Dict[str, list] = {}

    if project_initial:
        projection_result = weighted_projection(
            u_star=u,
            grid=grid,
            projection_data=projection_data,
            cfg=cfg,
        )
        u = projection_result.u_projected
        p_current = projection_result.phi

    t = 0.0

    while t < T - 1.0e-15:
        dt = choose_viscous_timestep(
            u=u,
            mu=mu,
            grid=grid,
            cfg=cfg,
            T_remaining=T - t,
            CFL=CFL,
            dt_max=dt_max,
        )

        u_star = viscous_predictor_explicit(
            u=u,
            mu=mu,
            dt=dt,
            grid=grid,
        )

        projection_result = weighted_projection(
            u_star=u_star,
            grid=grid,
            projection_data=projection_data,
            cfg=cfg,
        )

        u_new = projection_result.u_projected
        u_new[~grid.mask] = 0.0

        # Projection potential as pressure-like variable.
        p_current = projection_result.phi / dt if dt > 0.0 else projection_result.phi

        diag = compute_viscous_snapshot_diagnostics(
            u=u_new,
            projection_info=projection_result.info,
            mu=mu,
            dt=dt,
            time=t + dt,
            grid=grid,
        )

        _append_diag(diagnostics_history, diag)

        u = u_new
        t += dt

    diagnostics = _finalize_diag(diagnostics_history)

    return ViscousResult(
        u_final=u,
        p_final=p_current,
        diagnostics=diagnostics,
    )


# =============================================================================
# Printing
# =============================================================================

def print_viscous_summary(
    result: ViscousResult,
) -> None:
    """
    Print a compact summary of a viscous simulation.

    Parameters
    ----------
    result:
        ViscousResult object.
    """
    d = result.diagnostics

    print("=" * 70)
    print("Viscous solver summary")
    print("=" * 70)

    if "time" in d and len(d["time"]) > 0:
        print(f"final time                         = {d['time'][-1]:.8e}")
        print(f"number of time steps               = {len(d['time'])}")

    if "mu" in d and len(d["mu"]) > 0:
        print(f"mu                                 = {d['mu'][-1]:.8e}")

    if "weighted_energy" in d and len(d["weighted_energy"]) > 0:
        print(f"initial stored energy              = {d['weighted_energy'][0]:.8e}")
        print(f"final weighted energy              = {d['weighted_energy'][-1]:.8e}")

    if "divbu_l2" in d and len(d["divbu_l2"]) > 0:
        print(f"final ||div(bu)||_L2               = {d['divbu_l2'][-1]:.8e}")

    if "relative_divbu" in d and len(d["relative_divbu"]) > 0:
        print(f"final relative div(bu)             = {d['relative_divbu'][-1]:.8e}")

    if "projection_defect_after" in d and len(d["projection_defect_after"]) > 0:
        print(f"final projection defect after      = {d['projection_defect_after'][-1]:.8e}")

    if "projection_solver_residual" in d and len(d["projection_solver_residual"]) > 0:
        print(f"final projection solver residual   = {d['projection_solver_residual'][-1]:.8e}")

    if "projection_cg_iterations" in d and len(d["projection_cg_iterations"]) > 0:
        print(f"average projection CG iterations   = {np.mean(d['projection_cg_iterations']):.3f}")

    if "cfl_ratio" in d and len(d["cfl_ratio"]) > 0:
        print(f"max CFL ratio                      = {np.max(d['cfl_ratio']):.8e}")

    if "boundary_no_penetration_max" in d and len(d["boundary_no_penetration_max"]) > 0:
        print(f"boundary no-penetration max        = {d['boundary_no_penetration_max'][-1]:.8e}")

    print("=" * 70)