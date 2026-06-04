"""
Quick end-to-end test for the lake-equation vanishing-viscosity project.

This script checks the complete modular pipeline:

1. Configuration.
2. Grid and bathymetry.
3. Initial fields.
4. Singular elliptic reconstruction.
5. Inviscid solver.
6. Degenerate weighted projection.
7. Viscous solver.
8. Weighted velocity comparison.

Run from the project root directory:

    python scripts/run_quick_test.py

This is not a research-grade convergence study.
It is only a compact consistency test for the modular project.
"""

from pathlib import Path
import sys

import numpy as np


# =============================================================================
# Make project root importable
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =============================================================================
# Project imports
# =============================================================================

from lake.config import default_config, print_config_summary
from lake.grid import create_grid, print_grid_summary
from lake.operators import build_initial_fields, print_initial_diagnostics

from lake.elliptic import (
    assemble_singular_matrix,
    print_singular_matrix_summary,
    reconstruct_from_vorticity,
    reconstruction_diagnostics,
    print_reconstruction_diagnostics,
    reconstruction_sign_warning,
)

from lake.inviscid import (
    run_inviscid_solver,
    print_inviscid_summary,
)

from lake.projection import (
    build_projection_data,
    print_projection_matrix_summary,
    weighted_projection,
    print_projection_diagnostics,
)

from lake.viscous import (
    run_viscous_solver,
    print_viscous_summary,
)

from lake.diagnostics import (
    weighted_l2b_velocity_error,
    relative_weighted_l2b_velocity_error,
    build_vanishing_viscosity_table,
    print_vanishing_viscosity_table,
)


# =============================================================================
# Main quick test
# =============================================================================

def main() -> None:
    """
    Run a compact end-to-end test of the full modular project.
    """
    print()
    print("=" * 80)
    print("QUICK TEST — lake-equation vanishing-viscosity prototype")
    print("=" * 80)

    # -------------------------------------------------------------------------
    # 1. Configuration
    # -------------------------------------------------------------------------

    cfg = default_config()

    # Runtime-friendly defaults for quick testing.
    # If this is too slow, reduce cfg.N to 64 and T_limit to 0.02.
    cfg.N = 96
    cfg.T_inviscid = 0.05
    cfg.T_viscous = 0.05
    cfg.T_limit = 0.05
    cfg.mu_values = [1.0e-2, 5.0e-3]
    cfg.mu_plot_values = [1.0e-2]
    cfg.verbose = True

    print_config_summary(cfg)

    # -------------------------------------------------------------------------
    # 2. Grid and bathymetry
    # -------------------------------------------------------------------------

    grid = create_grid(cfg)
    print_grid_summary(grid, cfg)

    # -------------------------------------------------------------------------
    # 3. Initial fields
    # -------------------------------------------------------------------------

    fields = build_initial_fields(grid)
    print_initial_diagnostics(fields)

    # -------------------------------------------------------------------------
    # 4. Singular elliptic reconstruction
    # -------------------------------------------------------------------------

    ell_data = assemble_singular_matrix(grid)
    print_singular_matrix_summary(ell_data)

    rec = reconstruct_from_vorticity(
        omega=fields.omega0,
        grid=grid,
        elliptic_data=ell_data,
        cfg=cfg,
    )

    rec_diag = reconstruction_diagnostics(
        psi_reference=fields.psi0,
        u_reference=fields.u0,
        omega=fields.omega0,
        reconstruction=rec,
        grid=grid,
    )

    print_reconstruction_diagnostics(rec_diag)

    if reconstruction_sign_warning(fields.psi0, rec.psi, grid):
        print("WARNING: possible stream-function sign convention issue.")
    else:
        print("Sign convention check passed.")

    # -------------------------------------------------------------------------
    # 5. Inviscid solver
    # -------------------------------------------------------------------------

    result_inv = run_inviscid_solver(
        omega_initial=fields.omega0,
        grid=grid,
        elliptic_data=ell_data,
        cfg=cfg,
        T=cfg.T_limit,
        CFL=cfg.CFL_inviscid,
        dt_max=cfg.dt_max_inviscid,
    )

    print_inviscid_summary(result_inv)

    # -------------------------------------------------------------------------
    # 6. Degenerate weighted projection
    # -------------------------------------------------------------------------

    proj_data = build_projection_data(grid, cfg)
    print_projection_matrix_summary(proj_data)

    # Projection test on u0
    proj_u0 = weighted_projection(
        u_star=fields.u0,
        grid=grid,
        projection_data=proj_data,
        cfg=cfg,
    )

    print()
    print("Projection test on u0")
    print("-" * 80)
    print_projection_diagnostics(proj_u0)

    # Projection test on perturbed velocity
    rng = np.random.default_rng(123)

    perturbation = 0.05 * rng.standard_normal(size=fields.u0.shape)
    perturbation[~grid.mask] = 0.0

    u_star_perturbed = fields.u0 + perturbation
    u_star_perturbed[~grid.mask] = 0.0

    proj_pert = weighted_projection(
        u_star=u_star_perturbed,
        grid=grid,
        projection_data=proj_data,
        cfg=cfg,
    )

    print()
    print("Projection test on perturbed velocity")
    print("-" * 80)
    print_projection_diagnostics(proj_pert)

    # -------------------------------------------------------------------------
    # 7. Viscous solver for selected mu values
    # -------------------------------------------------------------------------

    viscous_runs = {}

    for mu in cfg.mu_values:
        print()
        print("-" * 80)
        print(f"Running viscous solver for mu = {mu:.4e}")
        print("-" * 80)

        result_visc = run_viscous_solver(
            u_initial=fields.u0,
            mu=float(mu),
            grid=grid,
            projection_data=proj_data,
            cfg=cfg,
            T=cfg.T_limit,
            CFL=cfg.CFL_viscous,
            dt_max=cfg.dt_max_viscous,
            project_initial=True,
        )

        print_viscous_summary(result_visc)

        viscous_runs[float(mu)] = result_visc

    # -------------------------------------------------------------------------
    # 8. Vanishing-viscosity velocity comparison
    # -------------------------------------------------------------------------

    u_mu_values = {
        float(mu): result.u_final
        for mu, result in viscous_runs.items()
    }

    viscous_diagnostics = {
        float(mu): result.diagnostics
        for mu, result in viscous_runs.items()
    }

    table = build_vanishing_viscosity_table(
        mu_values=cfg.mu_values,
        u_mu_values=u_mu_values,
        viscous_diagnostics=viscous_diagnostics,
        u_ref=result_inv.u_final,
        grid=grid,
    )

    print_vanishing_viscosity_table(table)

    # Direct explicit errors
    print()
    print("Direct weighted velocity errors")
    print("-" * 80)

    for mu, result in viscous_runs.items():
        err_abs = weighted_l2b_velocity_error(
            u_mu=result.u_final,
            u_ref=result_inv.u_final,
            grid=grid,
        )

        err_rel = relative_weighted_l2b_velocity_error(
            u_mu=result.u_final,
            u_ref=result_inv.u_final,
            grid=grid,
        )

        print(f"mu = {mu:.4e}")
        print(f"  ||u_mu - u||_L2_b       = {err_abs:.8e}")
        print(f"  relative weighted error = {err_rel:.8e}")

    # -------------------------------------------------------------------------
    # Final message
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("QUICK TEST COMPLETED")
    print("=" * 80)
    print("The modular project ran end-to-end.")
    print()
    print("Interpretation note:")
    print("  The observed vanishing-viscosity errors are qualitative.")
    print("  They combine viscosity error, h-error, dt-error, solver error,")
    print("  projection consistency error, and boundary/mask geometry error.")
    print("=" * 80)


if __name__ == "__main__":
    main()