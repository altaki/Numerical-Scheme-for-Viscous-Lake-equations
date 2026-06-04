"""
Configuration module for the lake-equation vanishing-viscosity project.

This file centralizes all physical, numerical, and solver parameters used by
the project.

The project studies:

1. Inviscid lake equations in vorticity-stream formulation.
2. Viscous lake equations in velocity formulation.
3. The vanishing-viscosity limit mu -> 0 in the weighted L^2_b norm.

Important mathematical conventions used throughout the project:

    grad_perp(psi) = (-partial_y psi, partial_x psi)

    curl(u) = partial_x u_y - partial_y u_x

    u = (1/b) grad_perp(psi)

    omega = curl(u) / b

    div((1/b) grad psi) = b omega

These conventions must not be changed in other modules.
"""

from dataclasses import dataclass, field
from typing import Tuple, List


@dataclass
class LakeConfig:
    """
    Main configuration class.

    This dataclass contains the numerical and physical parameters used by the
    whole project.

    The default values are chosen to keep the prototype reasonably fast.
    For more accurate experiments, increase N and refine dt, but expect longer
    runtimes.
    """

    # -------------------------------------------------------------------------
    # Grid and domain parameters
    # -------------------------------------------------------------------------

    #: Number of grid points in each direction on [-1,1]^2.
    #: N=96 is a robust default for a lightweight prototype.
    #: N=128 gives better resolution but may be slower.
    N: int = 96

    #: Domain bounds for the square containing the unit disk.
    x_min: float = -1.0
    x_max: float = 1.0
    y_min: float = -1.0
    y_max: float = 1.0

    #: Radius of the disk domain.
    radius: float = 1.0

    # -------------------------------------------------------------------------
    # Bathymetry parameters
    # -------------------------------------------------------------------------

    #: Bathymetry exponent:
    #:
    #:     b(x,y) = (1-r^2)^alpha
    #:
    #: inside the disk.
    #:
    #: Since near r=1,
    #:
    #:     1-r^2 ~ 2 dist(x, boundary),
    #:
    #: this behaves like dist(x,boundary)^alpha.
    #:
    #: The theoretical regime considered here is 0 < alpha < 1/2.
    alpha: float = 0.4

    #: Defensive floor used only for numerical division by b.
    #: This is NOT a structural regularization of the singular inviscid operator.
    eps_safe: float = 1.0e-14

    # -------------------------------------------------------------------------
    # Initial stream function parameters
    # -------------------------------------------------------------------------

    #: Center of the compactly supported initial stream-function bump.
    bump_center: Tuple[float, float] = (0.25, 0.0)

    #: Radius of the compactly supported bump.
    bump_radius: float = 0.25

    #: Threshold used for support-to-boundary diagnostics.
    support_threshold: float = 1.0e-3

    # -------------------------------------------------------------------------
    # Inviscid solver parameters
    # -------------------------------------------------------------------------

    #: Final time for a short inviscid demonstration.
    T_inviscid: float = 0.05

    #: CFL-like safety factor for semi-Lagrangian inviscid transport.
    CFL_inviscid: float = 0.5

    #: Maximum timestep for the inviscid solver.
    dt_max_inviscid: float = 2.0e-3

    # -------------------------------------------------------------------------
    # Viscous solver parameters
    # -------------------------------------------------------------------------

    #: Final time for short viscous demonstrations.
    T_viscous: float = 0.05

    #: CFL-like safety factor for the viscous velocity solver.
    CFL_viscous: float = 0.2

    #: Maximum timestep for the viscous solver.
    dt_max_viscous: float = 5.0e-4

    #: Viscosities used for quick viscous demonstrations.
    mu_plot_values: List[float] = field(
        default_factory=lambda: [1.0e-2, 5.0e-3]
    )

    # -------------------------------------------------------------------------
    # Vanishing-viscosity experiment parameters
    # -------------------------------------------------------------------------

    #: Common final time for comparing viscous and inviscid velocities.
    T_limit: float = 0.05

    #: Viscosity values used in the vanishing-viscosity experiment.
    #: Keep this list short for the prototype.
    mu_values: List[float] = field(
        default_factory=lambda: [1.0e-2, 5.0e-3, 2.0e-3]
    )

    #: Formal beta values used for plotting theoretical reference slopes:
    #:
    #:     mu^((1-beta)/2)
    #:
    #: These are only references unless the eta_mu boundary scaling and all
    #: numerical error budgets are controlled.
    beta_reference_values: List[float] = field(
        default_factory=lambda: [0.0, 0.5]
    )

    # -------------------------------------------------------------------------
    # Projection parameters
    # -------------------------------------------------------------------------

    #: If True, use epsilon_proj = h^2 by default in the degenerate projection.
    #: This stabilizes the projection but makes it approximate.
    use_projection_regularization: bool = True

    #: Multiplicative factor for epsilon_proj = factor * h^2.
    epsilon_proj_factor: float = 1.0

    # -------------------------------------------------------------------------
    # Sparse solver tolerances
    # -------------------------------------------------------------------------

    #: Relative tolerance for CG elliptic solves.
    cg_rtol: float = 1.0e-10

    #: Absolute tolerance for CG elliptic solves.
    cg_atol: float = 1.0e-12

    #: Maximum number of CG iterations.
    cg_maxiter: int = 5000

    # -------------------------------------------------------------------------
    # Plotting parameters
    # -------------------------------------------------------------------------

    #: Number of arrows along one direction for quiver plots.
    quiver_subsample_target: int = 20

    #: Default figure size for field plots.
    field_figsize: Tuple[float, float] = (6.0, 5.0)

    #: Default figure size for multiple-panel plots.
    panel_figsize: Tuple[float, float] = (12.0, 10.0)

    # -------------------------------------------------------------------------
    # Runtime / verbosity
    # -------------------------------------------------------------------------

    #: If True, print diagnostic summaries during simulations.
    verbose: bool = True

    #: If True, run extra consistency checks.
    run_extra_checks: bool = True


def default_config() -> LakeConfig:
    """
    Return the default project configuration.

    This helper is useful in notebooks:

        from lake.config import default_config
        cfg = default_config()
    """
    return LakeConfig()


def print_config_summary(cfg: LakeConfig) -> None:
    """
    Print a compact summary of the configuration.

    Parameters
    ----------
    cfg:
        LakeConfig instance.
    """
    print("=" * 70)
    print("Lake-equation vanishing-viscosity configuration")
    print("=" * 70)

    print("\nGrid/domain")
    print(f"  N                 = {cfg.N}")
    print(f"  domain square      = [{cfg.x_min}, {cfg.x_max}] x [{cfg.y_min}, {cfg.y_max}]")
    print(f"  disk radius        = {cfg.radius}")

    print("\nBathymetry")
    print(f"  alpha             = {cfg.alpha}")
    print(f"  eps_safe          = {cfg.eps_safe}")

    print("\nInitial data")
    print(f"  bump_center       = {cfg.bump_center}")
    print(f"  bump_radius       = {cfg.bump_radius}")
    print(f"  support_threshold = {cfg.support_threshold}")

    print("\nInviscid solver")
    print(f"  T_inviscid        = {cfg.T_inviscid}")
    print(f"  CFL_inviscid      = {cfg.CFL_inviscid}")
    print(f"  dt_max_inviscid   = {cfg.dt_max_inviscid}")

    print("\nViscous solver")
    print(f"  T_viscous         = {cfg.T_viscous}")
    print(f"  CFL_viscous       = {cfg.CFL_viscous}")
    print(f"  dt_max_viscous    = {cfg.dt_max_viscous}")
    print(f"  mu_plot_values    = {cfg.mu_plot_values}")

    print("\nVanishing-viscosity experiment")
    print(f"  T_limit           = {cfg.T_limit}")
    print(f"  mu_values         = {cfg.mu_values}")
    print(f"  beta refs         = {cfg.beta_reference_values}")

    print("\nProjection")
    print(f"  use regularization = {cfg.use_projection_regularization}")
    print(f"  eps factor         = {cfg.epsilon_proj_factor}")

    print("\nCG solver")
    print(f"  cg_rtol           = {cfg.cg_rtol}")
    print(f"  cg_atol           = {cfg.cg_atol}")
    print(f"  cg_maxiter        = {cfg.cg_maxiter}")

    print("=" * 70)