"""
Diagnostics for the lake-equation vanishing-viscosity project.

This module contains tools for:

- weighted velocity errors,
- relative weighted velocity errors,
- velocity diagnostics,
- log-log slope fitting,
- safe extraction of diagnostics from solver histories,
- vanishing-viscosity summary tables.

The main comparison is velocity-only:

    ||u_mu - u||_{L^2_b}
    =
    (int b |u_mu-u|^2 dx)^{1/2}.

This is the comparison used for the numerical vanishing-viscosity study.
"""

from typing import Dict, Any, Iterable, Tuple, List

import numpy as np

from .grid import LakeGrid
from .operators import (
    weighted_l2b_velocity_norm,
    weighted_kinetic_energy,
    ordinary_l2_vector_norm,
    masked_l2_norm,
    div_bu,
)


# =============================================================================
# Velocity errors
# =============================================================================

def weighted_l2b_velocity_error(
    u_mu: np.ndarray,
    u_ref: np.ndarray,
    grid: LakeGrid,
) -> float:
    """
    Compute the weighted velocity error:

        ||u_mu - u_ref||_{L^2_b}
        =
        sqrt(sum_mask b |u_mu-u_ref|^2 dx dy).

    Parameters
    ----------
    u_mu:
        Viscous velocity field.

    u_ref:
        Reference inviscid velocity field.

    grid:
        LakeGrid instance.

    Returns
    -------
    error:
        Weighted L2_b velocity error.
    """
    diff = u_mu - u_ref

    return weighted_l2b_velocity_norm(
        diff,
        grid.b,
        grid.mask,
        grid.dx,
        grid.dy,
    )


def relative_weighted_l2b_velocity_error(
    u_mu: np.ndarray,
    u_ref: np.ndarray,
    grid: LakeGrid,
) -> float:
    """
    Compute the relative weighted velocity error:

        ||u_mu-u_ref||_{L^2_b} / max(||u_ref||_{L^2_b}, 1e-14).

    Parameters
    ----------
    u_mu:
        Viscous velocity field.

    u_ref:
        Reference inviscid velocity field.

    grid:
        LakeGrid instance.

    Returns
    -------
    relative_error:
        Relative weighted L2_b velocity error.
    """
    error = weighted_l2b_velocity_error(
        u_mu,
        u_ref,
        grid,
    )

    ref_norm = weighted_l2b_velocity_norm(
        u_ref,
        grid.b,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    return error / max(ref_norm, 1.0e-14)


# =============================================================================
# Velocity diagnostics
# =============================================================================

def velocity_max_speed(
    u: np.ndarray,
    grid: LakeGrid,
) -> float:
    """
    Compute max |u| over the active mask.

    Parameters
    ----------
    u:
        Velocity field.

    grid:
        LakeGrid instance.

    Returns
    -------
    max_speed:
        Maximum velocity magnitude over the mask.
    """
    speed = np.sqrt(u[..., 0] ** 2 + u[..., 1] ** 2)

    if not np.any(grid.mask):
        return np.nan

    return float(np.max(speed[grid.mask]))


def divbu_l2_norm(
    u: np.ndarray,
    grid: LakeGrid,
) -> float:
    """
    Compute ||div(bu)||_L2.

    Parameters
    ----------
    u:
        Velocity field.

    grid:
        LakeGrid instance.

    Returns
    -------
    norm:
        L2 norm of weighted divergence.
    """
    div_val = div_bu(
        u,
        grid.b,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    return masked_l2_norm(
        div_val,
        grid.mask,
        grid.dx,
        grid.dy,
    )


def relative_divbu_norm(
    u: np.ndarray,
    grid: LakeGrid,
) -> float:
    """
    Compute a dimensionally scaled relative weighted-divergence defect:

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
        Relative weighted divergence defect.
    """
    div_l2 = divbu_l2_norm(
        u,
        grid,
    )

    bu_l2 = ordinary_l2_vector_norm(
        grid.b[..., None] * u,
        grid.mask,
        grid.dx,
        grid.dy,
    )

    return div_l2 / max(bu_l2 / grid.h, 1.0e-14)


def compute_velocity_diagnostics(
    u: np.ndarray,
    grid: LakeGrid,
) -> Dict[str, float]:
    """
    Compute standard velocity diagnostics.

    Parameters
    ----------
    u:
        Velocity field.

    grid:
        LakeGrid instance.

    Returns
    -------
    diagnostics:
        Dictionary containing:
        - weighted_energy
        - l2b_norm
        - divbu_l2
        - relative_divbu
        - max_speed
    """
    return {
        "weighted_energy": weighted_kinetic_energy(
            u,
            grid.b,
            grid.mask,
            grid.dx,
            grid.dy,
        ),
        "l2b_norm": weighted_l2b_velocity_norm(
            u,
            grid.b,
            grid.mask,
            grid.dx,
            grid.dy,
        ),
        "divbu_l2": divbu_l2_norm(u, grid),
        "relative_divbu": relative_divbu_norm(u, grid),
        "max_speed": velocity_max_speed(u, grid),
    }


# =============================================================================
# Log-log slope fitting
# =============================================================================

def fit_loglog_slope(
    x_values: Iterable[float],
    y_values: Iterable[float],
) -> Tuple[float, float]:
    """
    Fit a log-log slope:

        log(y) = slope log(x) + intercept.

    Nonpositive or nonfinite entries are ignored.

    Parameters
    ----------
    x_values:
        Positive x-values, e.g. mu values.

    y_values:
        Positive y-values, e.g. errors.

    Returns
    -------
    slope, intercept:
        Fit parameters.

        If fewer than two valid points are available, returns (nan, nan).
    """
    x = np.asarray(list(x_values), dtype=float)
    y = np.asarray(list(y_values), dtype=float)

    valid = (
        np.isfinite(x)
        & np.isfinite(y)
        & (x > 0.0)
        & (y > 0.0)
    )

    if np.sum(valid) < 2:
        return np.nan, np.nan

    coeffs = np.polyfit(
        np.log(x[valid]),
        np.log(y[valid]),
        1,
    )

    return float(coeffs[0]), float(coeffs[1])


def evaluate_power_law(
    x: np.ndarray,
    slope: float,
    intercept: float,
) -> np.ndarray:
    """
    Evaluate exp(intercept) * x^slope.

    Parameters
    ----------
    x:
        Positive x-values.

    slope:
        Power-law slope.

    intercept:
        Log-space intercept.

    Returns
    -------
    y:
        Power-law values.
    """
    return np.exp(intercept) * x**slope


# =============================================================================
# Safe diagnostic extraction
# =============================================================================

def get_last_diag_value(
    diagnostics: Dict[str, Any],
    key: str,
    default: float = np.nan,
) -> float:
    """
    Safely get the last value of a diagnostic array.

    Parameters
    ----------
    diagnostics:
        Dictionary of diagnostic arrays.

    key:
        Key to extract.

    default:
        Value returned if key is absent or empty.

    Returns
    -------
    value:
        Last diagnostic value or default.
    """
    if key not in diagnostics:
        return default

    value = diagnostics[key]

    try:
        arr = np.asarray(value)
        if arr.size == 0:
            return default
        return float(arr[-1])
    except Exception:
        return default


def get_mean_diag_value(
    diagnostics: Dict[str, Any],
    key: str,
    default: float = np.nan,
) -> float:
    """
    Safely compute the mean value of a diagnostic array.

    Parameters
    ----------
    diagnostics:
        Dictionary of diagnostic arrays.

    key:
        Key to extract.

    default:
        Value returned if key is absent or empty.

    Returns
    -------
    mean_value:
        Mean diagnostic value or default.
    """
    if key not in diagnostics:
        return default

    value = diagnostics[key]

    try:
        arr = np.asarray(value)
        if arr.size == 0:
            return default
        return float(np.mean(arr))
    except Exception:
        return default


def get_max_diag_value(
    diagnostics: Dict[str, Any],
    key: str,
    default: float = np.nan,
) -> float:
    """
    Safely compute the maximum value of a diagnostic array.

    Parameters
    ----------
    diagnostics:
        Dictionary of diagnostic arrays.

    key:
        Key to extract.

    default:
        Value returned if key is absent or empty.

    Returns
    -------
    max_value:
        Maximum diagnostic value or default.
    """
    if key not in diagnostics:
        return default

    value = diagnostics[key]

    try:
        arr = np.asarray(value)
        if arr.size == 0:
            return default
        return float(np.max(arr))
    except Exception:
        return default


# =============================================================================
# Vanishing-viscosity summaries
# =============================================================================

def summarize_viscous_result(
    u_mu: np.ndarray,
    diagnostics: Dict[str, Any],
    grid: LakeGrid,
) -> Dict[str, float]:
    """
    Summarize one viscous result.

    Parameters
    ----------
    u_mu:
        Final viscous velocity.

    diagnostics:
        Viscous solver diagnostics.

    grid:
        LakeGrid instance.

    Returns
    -------
    summary:
        Dictionary of final and aggregate diagnostics.
    """
    vel_diag = compute_velocity_diagnostics(
        u_mu,
        grid,
    )

    summary = {
        "weighted_energy": vel_diag["weighted_energy"],
        "l2b_norm": vel_diag["l2b_norm"],
        "divbu_l2": vel_diag["divbu_l2"],
        "relative_divbu": vel_diag["relative_divbu"],
        "max_speed": vel_diag["max_speed"],
        "projection_defect_after": get_last_diag_value(
            diagnostics,
            "projection_defect_after",
        ),
        "projection_solver_residual": get_last_diag_value(
            diagnostics,
            "projection_solver_residual",
        ),
        "average_projection_cg_iterations": get_mean_diag_value(
            diagnostics,
            "projection_cg_iterations",
        ),
        "max_cfl_ratio": get_max_diag_value(
            diagnostics,
            "cfl_ratio",
        ),
        "boundary_no_penetration_max": get_last_diag_value(
            diagnostics,
            "boundary_no_penetration_max",
        ),
    }

    return summary


def build_vanishing_viscosity_table(
    mu_values: Iterable[float],
    u_mu_values: Dict[float, np.ndarray],
    viscous_diagnostics: Dict[float, Dict[str, Any]],
    u_ref: np.ndarray,
    grid: LakeGrid,
) -> Dict[str, np.ndarray]:
    """
    Build arrays for a vanishing-viscosity comparison.

    Parameters
    ----------
    mu_values:
        Iterable of viscosity values.

    u_mu_values:
        Dictionary mapping mu -> final viscous velocity.

    viscous_diagnostics:
        Dictionary mapping mu -> viscous diagnostic dictionary.

    u_ref:
        Reference inviscid velocity.

    grid:
        LakeGrid instance.

    Returns
    -------
    table:
        Dictionary of numpy arrays containing:
        - mu
        - error_abs
        - error_rel
        - divbu_l2
        - relative_divbu
        - projection_defect_after
        - projection_solver_residual
        - average_projection_cg_iterations
        - max_cfl_ratio
        - weighted_energy
    """
    mu_list: List[float] = []
    err_abs_list: List[float] = []
    err_rel_list: List[float] = []
    div_list: List[float] = []
    rel_div_list: List[float] = []
    proj_def_list: List[float] = []
    proj_res_list: List[float] = []
    cg_iter_list: List[float] = []
    cfl_list: List[float] = []
    energy_list: List[float] = []

    for mu in mu_values:
        mu_float = float(mu)

        if mu_float not in u_mu_values:
            continue

        u_mu = u_mu_values[mu_float]
        diag = viscous_diagnostics.get(mu_float, {})

        summary = summarize_viscous_result(
            u_mu,
            diag,
            grid,
        )

        err_abs = weighted_l2b_velocity_error(
            u_mu,
            u_ref,
            grid,
        )

        err_rel = relative_weighted_l2b_velocity_error(
            u_mu,
            u_ref,
            grid,
        )

        mu_list.append(mu_float)
        err_abs_list.append(err_abs)
        err_rel_list.append(err_rel)
        div_list.append(summary["divbu_l2"])
        rel_div_list.append(summary["relative_divbu"])
        proj_def_list.append(summary["projection_defect_after"])
        proj_res_list.append(summary["projection_solver_residual"])
        cg_iter_list.append(summary["average_projection_cg_iterations"])
        cfl_list.append(summary["max_cfl_ratio"])
        energy_list.append(summary["weighted_energy"])

    table = {
        "mu": np.asarray(mu_list, dtype=float),
        "error_abs": np.asarray(err_abs_list, dtype=float),
        "error_rel": np.asarray(err_rel_list, dtype=float),
        "divbu_l2": np.asarray(div_list, dtype=float),
        "relative_divbu": np.asarray(rel_div_list, dtype=float),
        "projection_defect_after": np.asarray(proj_def_list, dtype=float),
        "projection_solver_residual": np.asarray(proj_res_list, dtype=float),
        "average_projection_cg_iterations": np.asarray(cg_iter_list, dtype=float),
        "max_cfl_ratio": np.asarray(cfl_list, dtype=float),
        "weighted_energy": np.asarray(energy_list, dtype=float),
    }

    slope_abs, intercept_abs = fit_loglog_slope(
        table["mu"],
        table["error_abs"],
    )

    slope_rel, intercept_rel = fit_loglog_slope(
        table["mu"],
        table["error_rel"],
    )

    table["slope_abs"] = np.asarray([slope_abs])
    table["intercept_abs"] = np.asarray([intercept_abs])
    table["slope_rel"] = np.asarray([slope_rel])
    table["intercept_rel"] = np.asarray([intercept_rel])

    return table


def print_vanishing_viscosity_table(
    table: Dict[str, np.ndarray],
) -> None:
    """
    Print a vanishing-viscosity summary table.

    Parameters
    ----------
    table:
        Dictionary returned by build_vanishing_viscosity_table.
    """
    mu = table.get("mu", np.asarray([]))

    print("=" * 110)
    print("Vanishing-viscosity velocity comparison")
    print("=" * 110)

    if mu.size == 0:
        print("No data available.")
        print("=" * 110)
        return

    header = (
        f"{'mu':>12} | "
        f"{'err_abs':>14} | "
        f"{'err_rel':>14} | "
        f"{'divbu':>14} | "
        f"{'proj_def':>14} | "
        f"{'solver_res':>14} | "
        f"{'max_CFL':>10}"
    )

    print(header)
    print("-" * len(header))

    for i in range(mu.size):
        print(
            f"{table['mu'][i]:12.4e} | "
            f"{table['error_abs'][i]:14.6e} | "
            f"{table['error_rel'][i]:14.6e} | "
            f"{table['divbu_l2'][i]:14.6e} | "
            f"{table['projection_defect_after'][i]:14.6e} | "
            f"{table['projection_solver_residual'][i]:14.6e} | "
            f"{table['max_cfl_ratio'][i]:10.4e}"
        )

    slope_abs = table.get("slope_abs", np.asarray([np.nan]))[0]
    slope_rel = table.get("slope_rel", np.asarray([np.nan]))[0]

    print("-" * len(header))
    print(f"fitted absolute-error slope = {slope_abs:.8e}")
    print(f"fitted relative-error slope = {slope_rel:.8e}")
    print()
    print("Reference theoretical slopes:")
    print("  beta = 0.0 -> (1-beta)/2 = 0.5")
    print("  beta = 0.5 -> (1-beta)/2 = 0.25")
    print()
    print("Interpret slopes qualitatively unless h, dt, solver, projection,")
    print("boundary, and eta_mu-scaling errors are controlled.")
    print("=" * 110)
