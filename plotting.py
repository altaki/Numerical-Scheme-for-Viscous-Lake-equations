"""
Plotting utilities for the lake-equation vanishing-viscosity project.

This module contains plotting functions for:

- bathymetry and initial data,
- inviscid vorticity-stream solutions,
- viscous velocity solutions,
- time diagnostics,
- vanishing-viscosity error curves.

The plotting functions are deliberately lightweight and use only matplotlib.

All field plots mask values outside the disk by replacing them with NaN.
"""

from typing import Dict, Any, Iterable, Optional, Sequence

import numpy as np
import matplotlib.pyplot as plt

from .grid import LakeGrid
from .operators import (
    div_bu,
    curl_vector,
)
from .diagnostics import evaluate_power_law


# =============================================================================
# Basic plot helpers
# =============================================================================

def mask_to_nan(
    arr: np.ndarray,
    grid: LakeGrid,
) -> np.ndarray:
    """
    Return arr inside the active mask and NaN outside.

    Parameters
    ----------
    arr:
        Scalar field.

    grid:
        LakeGrid instance.

    Returns
    -------
    masked_arr:
        Array with NaN outside the active disk.
    """
    return np.where(grid.mask, arr, np.nan)


def velocity_magnitude(
    u: np.ndarray,
) -> np.ndarray:
    """
    Compute velocity magnitude |u|.

    Parameters
    ----------
    u:
        Vector field with shape (..., 2).

    Returns
    -------
    mag:
        Scalar field |u|.
    """
    return np.sqrt(u[..., 0] ** 2 + u[..., 1] ** 2)


def add_colorbar(
    fig: plt.Figure,
    ax: plt.Axes,
    im,
) -> None:
    """
    Add a compact colorbar to an axis.
    """
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def plot_scalar_field(
    field: np.ndarray,
    grid: LakeGrid,
    title: str = "",
    cmap: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    colorbar: bool = True,
):
    """
    Plot a scalar field on the masked disk.

    Parameters
    ----------
    field:
        Scalar field.

    grid:
        LakeGrid instance.

    title:
        Plot title.

    cmap:
        Optional matplotlib colormap.

    ax:
        Optional axis. If None, creates a new figure.

    colorbar:
        If True, add a colorbar.

    Returns
    -------
    ax:
        Axis object.
    """
    created_fig = False

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
        created_fig = True
    else:
        fig = ax.figure

    im = ax.imshow(
        mask_to_nan(field, grid),
        origin="lower",
        extent=[grid.x.min(), grid.x.max(), grid.y.min(), grid.y.max()],
        cmap=cmap,
    )

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    if colorbar:
        add_colorbar(fig, ax, im)

    if created_fig:
        plt.show()

    return ax


def plot_quiver(
    u: np.ndarray,
    grid: LakeGrid,
    title: str = "velocity quiver",
    ax: Optional[plt.Axes] = None,
    target: int = 20,
):
    """
    Plot a coarse quiver plot of a velocity field.

    Parameters
    ----------
    u:
        Velocity field.

    grid:
        LakeGrid instance.

    title:
        Plot title.

    ax:
        Optional axis.

    target:
        Approximate number of arrows per direction.

    Returns
    -------
    ax:
        Axis object.
    """
    created_fig = False

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
        created_fig = True

    step = max(grid.x.size // target, 1)

    Xs = grid.X[::step, ::step]
    Ys = grid.Y[::step, ::step]
    Us = u[::step, ::step, 0]
    Vs = u[::step, ::step, 1]
    Ms = grid.mask[::step, ::step]

    ax.quiver(
        Xs[Ms],
        Ys[Ms],
        Us[Ms],
        Vs[Ms],
    )

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    if created_fig:
        plt.show()

    return ax


# =============================================================================
# Initial data plots
# =============================================================================

def plot_initial_data(
    grid: LakeGrid,
    omega0: np.ndarray,
    u0: np.ndarray,
) -> None:
    """
    Plot bathymetry, initial stream function, initial vorticity and velocity.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    omega0:
        Initial vorticity.

    u0:
        Initial velocity.
    """
    u0_mag = velocity_magnitude(u0)

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12, 10),
        constrained_layout=True,
    )

    plot_scalar_field(
        grid.b,
        grid,
        title="Bathymetry b",
        ax=axes[0, 0],
    )

    plot_scalar_field(
        grid.psi0,
        grid,
        title="Initial stream function psi0",
        ax=axes[0, 1],
    )

    plot_scalar_field(
        omega0,
        grid,
        title="Initial vorticity omega0",
        ax=axes[1, 0],
    )

    plot_scalar_field(
        u0_mag,
        grid,
        title="Initial velocity magnitude |u0|",
        ax=axes[1, 1],
    )

    plt.show()


# =============================================================================
# Inviscid plots
# =============================================================================

def plot_inviscid_solution(
    omega_initial: np.ndarray,
    omega_final: np.ndarray,
    psi_final: np.ndarray,
    u_final: np.ndarray,
    grid: LakeGrid,
) -> None:
    """
    Plot inviscid initial/final vorticity, stream function and velocity.

    Parameters
    ----------
    omega_initial:
        Initial vorticity.

    omega_final:
        Final vorticity.

    psi_final:
        Final stream function.

    u_final:
        Final velocity.
    """
    u_mag = velocity_magnitude(u_final)

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12, 10),
        constrained_layout=True,
    )

    plot_scalar_field(
        omega_initial,
        grid,
        title="Initial vorticity omega0",
        ax=axes[0, 0],
    )

    plot_scalar_field(
        omega_final,
        grid,
        title="Final inviscid vorticity omega",
        ax=axes[0, 1],
    )

    plot_scalar_field(
        psi_final,
        grid,
        title="Final stream function psi",
        ax=axes[1, 0],
    )

    plot_scalar_field(
        u_mag,
        grid,
        title="Final inviscid velocity magnitude |u|",
        ax=axes[1, 1],
    )

    plt.show()


def plot_inviscid_diagnostics(
    diagnostics: Dict[str, np.ndarray],
) -> None:
    """
    Plot time histories from the inviscid solver.

    Parameters
    ----------
    diagnostics:
        Diagnostics dictionary from run_inviscid_solver.
    """
    if "time" not in diagnostics:
        print("No time array found in inviscid diagnostics.")
        return

    t = diagnostics["time"]

    fig, axes = plt.subplots(
        3,
        2,
        figsize=(12, 12),
        constrained_layout=True,
    )

    if "min_omega" in diagnostics and "max_omega" in diagnostics:
        axes[0, 0].plot(t, diagnostics["min_omega"], label="min omega")
        axes[0, 0].plot(t, diagnostics["max_omega"], label="max omega")
        axes[0, 0].legend()
    axes[0, 0].set_title("min/max omega")
    axes[0, 0].set_xlabel("time")

    if "support_distance" in diagnostics:
        axes[0, 1].plot(t, diagnostics["support_distance"])
    axes[0, 1].set_title("support-to-boundary distance")
    axes[0, 1].set_xlabel("time")

    if "divbu_l2" in diagnostics:
        axes[1, 0].plot(t, diagnostics["divbu_l2"])
    axes[1, 0].set_title(r"$||div(bu)||_{L^2}$")
    axes[1, 0].set_xlabel("time")

    if "vorticity_residual_l2" in diagnostics:
        axes[1, 1].plot(t, diagnostics["vorticity_residual_l2"])
    axes[1, 1].set_title(r"$||curl(u)/b-\omega||_{L^2}$")
    axes[1, 1].set_xlabel("time")

    if "elliptic_residual" in diagnostics:
        axes[2, 0].plot(t, diagnostics["elliptic_residual"])
    axes[2, 0].set_title("elliptic residual")
    axes[2, 0].set_xlabel("time")

    if "cfl_ratio" in diagnostics:
        axes[2, 1].plot(t, diagnostics["cfl_ratio"])
    axes[2, 1].set_title("CFL ratio")
    axes[2, 1].set_xlabel("time")

    plt.show()


# =============================================================================
# Viscous plots
# =============================================================================

def plot_viscous_solution(
    u_mu: np.ndarray,
    p_mu: np.ndarray,
    grid: LakeGrid,
    title_suffix: str = "",
) -> None:
    """
    Plot viscous velocity magnitude, quiver, pressure/projection potential and div(bu).

    Parameters
    ----------
    u_mu:
        Viscous velocity.

    p_mu:
        Pressure/projection-potential approximation.

    grid:
        LakeGrid instance.

    title_suffix:
        Extra title text, e.g. "mu=1e-2".
    """
    speed = velocity_magnitude(u_mu)

    div_val = div_bu(
        u_mu,
        grid.b,
        grid.dx,
        grid.dy,
        grid.mask,
    )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12, 10),
        constrained_layout=True,
    )

    plot_scalar_field(
        speed,
        grid,
        title=f"|u_mu| {title_suffix}",
        ax=axes[0, 0],
    )

    plot_quiver(
        u_mu,
        grid,
        title=f"velocity quiver {title_suffix}",
        ax=axes[0, 1],
    )

    plot_scalar_field(
        p_mu,
        grid,
        title=f"pressure/projection potential {title_suffix}",
        ax=axes[1, 0],
    )

    plot_scalar_field(
        div_val,
        grid,
        title=f"div(bu_mu) {title_suffix}",
        ax=axes[1, 1],
    )

    plt.show()


def plot_viscous_diagnostics(
    diagnostics: Dict[str, np.ndarray],
    title_suffix: str = "",
) -> None:
    """
    Plot time histories from the viscous solver.

    Parameters
    ----------
    diagnostics:
        Diagnostics dictionary from run_viscous_solver.

    title_suffix:
        Optional text for titles.
    """
    if "time" not in diagnostics:
        print("No time array found in viscous diagnostics.")
        return

    t = diagnostics["time"]

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(8, 10),
        constrained_layout=True,
    )

    if "weighted_energy" in diagnostics:
        axes[0].plot(t, diagnostics["weighted_energy"])
    axes[0].set_title(f"Weighted kinetic energy {title_suffix}")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel("energy")

    if "projection_defect_after" in diagnostics:
        axes[1].plot(t, diagnostics["projection_defect_after"])
    axes[1].set_title(f"Projection defect after {title_suffix}")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("defect")

    if "projection_solver_residual" in diagnostics:
        axes[2].plot(t, diagnostics["projection_solver_residual"])
    axes[2].set_title(f"Projection solver residual {title_suffix}")
    axes[2].set_xlabel("time")
    axes[2].set_ylabel("residual")

    plt.show()


# =============================================================================
# Vanishing-viscosity plots
# =============================================================================

def _valid_positive(
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """
    Boolean mask for positive finite pairs.
    """
    return (
        np.isfinite(x)
        & np.isfinite(y)
        & (x > 0.0)
        & (y > 0.0)
    )


def plot_vanishing_errors(
    table: Dict[str, np.ndarray],
    beta_reference_values: Sequence[float] = (0.0, 0.5),
) -> None:
    """
    Plot weighted velocity errors versus viscosity.

    Parameters
    ----------
    table:
        Dictionary returned by build_vanishing_viscosity_table.

    beta_reference_values:
        Beta values for reference slopes:

            mu^((1-beta)/2).
    """
    mu = table.get("mu", np.asarray([]))
    err_abs = table.get("error_abs", np.asarray([]))
    err_rel = table.get("error_rel", np.asarray([]))

    slope_abs = table.get("slope_abs", np.asarray([np.nan]))[0]
    intercept_abs = table.get("intercept_abs", np.asarray([np.nan]))[0]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13, 5),
        constrained_layout=True,
    )

    # Absolute error
    valid_abs = _valid_positive(mu, err_abs)

    if np.any(valid_abs):
        mu_abs = mu[valid_abs]
        err_abs_valid = err_abs[valid_abs]

        axes[0].loglog(
            mu_abs,
            err_abs_valid,
            "o-",
            label=r"$||u_\mu-u||_{L^2_b}$",
        )

        if np.isfinite(slope_abs) and np.isfinite(intercept_abs):
            mu_fit = np.exp(
                np.linspace(
                    np.log(np.min(mu_abs)),
                    np.log(np.max(mu_abs)),
                    100,
                )
            )
            err_fit = evaluate_power_law(
                mu_fit,
                slope_abs,
                intercept_abs,
            )
            axes[0].loglog(
                mu_fit,
                err_fit,
                "--",
                label=f"fit slope={slope_abs:.3f}",
            )

        # Reference slopes through the first point.
        mu0 = mu_abs[0]
        e0 = err_abs_valid[0]

        for beta in beta_reference_values:
            ref_slope = (1.0 - beta) / 2.0
            ref_values = e0 * (mu_abs / mu0) ** ref_slope
            axes[0].loglog(
                mu_abs,
                ref_values,
                ":",
                label=rf"reference $\mu^{{{ref_slope:.2f}}}$",
            )

    axes[0].set_title("Absolute weighted velocity error")
    axes[0].set_xlabel(r"$\mu$")
    axes[0].set_ylabel(r"$||u_\mu-u||_{L^2_b}$")
    axes[0].grid(True, which="both")
    axes[0].legend()

    # Relative error
    valid_rel = _valid_positive(mu, err_rel)

    if np.any(valid_rel):
        axes[1].loglog(
            mu[valid_rel],
            err_rel[valid_rel],
            "o-",
            label="relative error",
        )

    axes[1].set_title("Relative weighted velocity error")
    axes[1].set_xlabel(r"$\mu$")
    axes[1].set_ylabel("relative error")
    axes[1].grid(True, which="both")
    axes[1].legend()

    plt.show()


def plot_error_budget(
    table: Dict[str, np.ndarray],
) -> None:
    """
    Plot error-budget indicators versus viscosity.

    Parameters
    ----------
    table:
        Dictionary returned by build_vanishing_viscosity_table.
    """
    mu = table.get("mu", np.asarray([]))

    indicators = [
        ("error_abs", "velocity error"),
        ("divbu_l2", "div defect"),
        ("projection_defect_after", "projection defect"),
        ("projection_solver_residual", "solver residual"),
    ]

    plt.figure(figsize=(8, 6))

    plotted = False

    for key, label in indicators:
        if key not in table:
            continue

        y = table[key]
        valid = _valid_positive(mu, y)

        if np.any(valid):
            plt.loglog(
                mu[valid],
                y[valid],
                "o-",
                label=label,
            )
            plotted = True

    if plotted:
        plt.xlabel(r"$\mu$")
        plt.ylabel("indicator")
        plt.title("Error-budget indicators")
        plt.grid(True, which="both")
        plt.legend()
        plt.show()
    else:
        plt.close()
        print("No positive finite error-budget indicators to plot.")


def plot_vanishing_field_comparison(
    u_ref: np.ndarray,
    u_mu: np.ndarray,
    grid: LakeGrid,
    mu: float,
) -> None:
    """
    Compare inviscid and viscous final velocity fields.

    Parameters
    ----------
    u_ref:
        Reference inviscid velocity.

    u_mu:
        Viscous velocity.

    grid:
        LakeGrid instance.

    mu:
        Viscosity value.
    """
    speed_ref = velocity_magnitude(u_ref)
    speed_mu = velocity_magnitude(u_mu)
    speed_diff = velocity_magnitude(u_mu - u_ref)

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(15, 5),
        constrained_layout=True,
    )

    plot_scalar_field(
        speed_ref,
        grid,
        title="|u inviscid|",
        ax=axes[0],
    )

    plot_scalar_field(
        speed_mu,
        grid,
        title=rf"$|u_\mu|$, $\mu={mu:.2e}$",
        ax=axes[1],
    )

    plot_scalar_field(
        speed_diff,
        grid,
        title=rf"$|u_\mu-u|$, $\mu={mu:.2e}$",
        ax=axes[2],
    )

    plt.show()


def plot_multiple_viscous_histories(
    viscous_results: Dict[float, Any],
    key: str,
    title: str,
    ylabel: str,
) -> None:
    """
    Plot a given diagnostic history for multiple viscosity values.

    Parameters
    ----------
    viscous_results:
        Dictionary mapping mu to ViscousResult.

    key:
        Diagnostic key to plot.

    title:
        Plot title.

    ylabel:
        y-axis label.
    """
    plt.figure(figsize=(8, 5))

    plotted = False

    for mu, result in viscous_results.items():
        diag = result.diagnostics

        if "time" not in diag or key not in diag:
            continue

        plt.plot(
            diag["time"],
            diag[key],
            label=rf"$\mu={mu:.2e}$",
        )
        plotted = True

    if plotted:
        plt.xlabel("time")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.legend()
        plt.show()
    else:
        plt.close()
        print(f"No diagnostic history found for key={key!r}.")