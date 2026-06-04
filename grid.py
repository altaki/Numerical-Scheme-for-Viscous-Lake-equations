"""
Grid and bathymetry construction for the lake-equation project.

This module defines the masked Cartesian approximation of the disk domain,
the degenerate bathymetry b, the defensive coefficient b_safe, and the
compactly supported initial stream function psi0.

Domain:

    Omega = { (x,y) : x^2 + y^2 < R^2 }

where R = cfg.radius.

Bathymetry:

    b(x,y) = (R^2-r^2)^alpha,   r = sqrt(x^2+y^2)

inside the disk, and b=0 outside.

Near r=R:

    R^2-r^2 = (R-r)(R+r) ~ 2R dist(x, boundary),

so b behaves like dist(x,boundary)^alpha.

The default radius is R=1 and the default alpha=0.4, which is in the
regime 0 < alpha < 1/2.
"""

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from .config import LakeConfig


@dataclass
class LakeGrid:
    """
    Container for the grid, mask, bathymetry, and initial stream function.

    Attributes
    ----------
    radius:
        Radius of the disk domain.

    x, y:
        One-dimensional Cartesian grids.

    X, Y:
        Two-dimensional meshgrid arrays.

    r:
        Radial distance sqrt(X^2+Y^2).

    dx, dy, h:
        Grid spacings.

    mask:
        Boolean mask for the active disk nodes.

    dist_to_boundary:
        Distance to the circular boundary inside the active mask.

    boundary_band_topological:
        Active nodes with at least one inactive 4-neighbor.

    boundary_band_1h, boundary_band_2h, boundary_band_3h:
        Geometric boundary bands of thickness h, 2h, and 3h.

    boundary_band:
        Default diagnostic boundary band, currently boundary_band_3h.

    interior_unknowns:
        Active nodes excluding the topological boundary band.

    b:
        Bathymetry.

    b_safe:
        Defensive version of b used only for division.

    psi0:
        Compactly supported initial stream function.
    """

    radius: float

    x: np.ndarray
    y: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    r: np.ndarray

    dx: float
    dy: float
    h: float

    mask: np.ndarray

    dist_to_boundary: np.ndarray
    boundary_band_topological: np.ndarray
    boundary_band_1h: np.ndarray
    boundary_band_2h: np.ndarray
    boundary_band_3h: np.ndarray
    boundary_band: np.ndarray

    interior_unknowns: np.ndarray

    b: np.ndarray
    b_safe: np.ndarray

    psi0: np.ndarray


def create_cartesian_grid(
    cfg: LakeConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    """
    Create the Cartesian grid and radial coordinate.

    Parameters
    ----------
    cfg:
        LakeConfig instance.

    Returns
    -------
    x, y, X, Y, r, dx, dy, h
    """
    x = np.linspace(cfg.x_min, cfg.x_max, cfg.N)
    y = np.linspace(cfg.y_min, cfg.y_max, cfg.N)

    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    h = max(dx, dy)

    X, Y = np.meshgrid(x, y, indexing="xy")
    r = np.sqrt(X**2 + Y**2)

    return x, y, X, Y, r, dx, dy, h


def create_disk_mask(
    r: np.ndarray,
    cfg: LakeConfig,
) -> np.ndarray:
    """
    Create the active mask for the disk.

    Parameters
    ----------
    r:
        Radial coordinate.

    cfg:
        LakeConfig instance.

    Returns
    -------
    mask:
        Boolean array for the disk r < radius.
    """
    return r < cfg.radius


def compute_topological_boundary_band(
    mask: np.ndarray,
) -> np.ndarray:
    """
    Compute the topological boundary band of a masked Cartesian disk.

    A node belongs to the topological boundary band if:

    - it is inside the active mask;
    - at least one of its four Cartesian neighbors is outside the mask.

    Parameters
    ----------
    mask:
        Boolean active mask.

    Returns
    -------
    boundary_band_topological:
        Boolean array marking boundary-adjacent active nodes.
    """
    mask_pad = np.pad(mask.astype(bool), 1, mode="constant", constant_values=False)

    up = mask_pad[:-2, 1:-1]
    down = mask_pad[2:, 1:-1]
    left = mask_pad[1:-1, :-2]
    right = mask_pad[1:-1, 2:]

    boundary_band_topological = mask & (~up | ~down | ~left | ~right)

    return boundary_band_topological


def compute_geometric_boundary_bands(
    r: np.ndarray,
    mask: np.ndarray,
    h: float,
    cfg: LakeConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute geometric boundary bands of thickness h, 2h, and 3h.

    For a disk of radius R:

        dist_to_boundary = R - r.

    Parameters
    ----------
    r:
        Radial coordinate.

    mask:
        Active disk mask.

    h:
        Grid spacing.

    cfg:
        LakeConfig instance.

    Returns
    -------
    dist_to_boundary, boundary_band_1h, boundary_band_2h, boundary_band_3h
    """
    dist_to_boundary = np.full_like(r, np.nan, dtype=float)
    dist_to_boundary[mask] = cfg.radius - r[mask]

    boundary_band_1h = mask & (dist_to_boundary <= 1.0 * h)
    boundary_band_2h = mask & (dist_to_boundary <= 2.0 * h)
    boundary_band_3h = mask & (dist_to_boundary <= 3.0 * h)

    return dist_to_boundary, boundary_band_1h, boundary_band_2h, boundary_band_3h


def compute_bathymetry(
    r: np.ndarray,
    mask: np.ndarray,
    cfg: LakeConfig,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute bathymetry b and defensive b_safe.

    For a disk of radius R:

        b = (R^2-r^2)^alpha inside the disk,
        b = 0 outside the disk.

    Near r=R:

        R^2-r^2 = (R-r)(R+r) ~ 2R dist(x,boundary),

    so b behaves like dist(x,boundary)^alpha.

    Defensive coefficient:

        b_safe = max(b, eps_safe) inside the mask,
        b_safe = 1 outside the mask.

    b_safe is only used to prevent numerical division by zero. It should not be
    interpreted as a structural regularization of the singular operator 1/b.

    Parameters
    ----------
    r:
        Radial coordinate.

    mask:
        Active disk mask.

    cfg:
        LakeConfig instance.

    Returns
    -------
    b, b_safe
    """
    b = np.zeros_like(r, dtype=float)

    b[mask] = np.maximum(
        cfg.radius**2 - r[mask] ** 2,
        0.0,
    ) ** cfg.alpha

    b[~mask] = 0.0

    b_safe = np.where(
        mask,
        np.maximum(b, cfg.eps_safe),
        1.0,
    )

    return b, b_safe


def create_initial_stream_function(
    X: np.ndarray,
    Y: np.ndarray,
    mask: np.ndarray,
    cfg: LakeConfig,
) -> np.ndarray:
    """
    Create the compactly supported initial stream function psi0.

    The stream function is a smooth bump supported away from the boundary:

        center = cfg.bump_center
        rho    = cfg.bump_radius

        s = ((X-center_x)^2 + (Y-center_y)^2)/rho^2

        psi0 = exp(-1/(1-s)) if s < 1,
        psi0 = 0 otherwise.

    Parameters
    ----------
    X, Y:
        Meshgrid arrays.

    mask:
        Active disk mask.

    cfg:
        LakeConfig instance.

    Returns
    -------
    psi0:
        Initial stream function.
    """
    cx, cy = cfg.bump_center
    rho = cfg.bump_radius

    s = ((X - cx) ** 2 + (Y - cy) ** 2) / rho**2

    psi0 = np.zeros_like(X, dtype=float)

    inside_bump = s < 1.0
    psi0[inside_bump] = np.exp(-1.0 / (1.0 - s[inside_bump]))

    psi0[~mask] = 0.0

    return psi0


def create_grid(
    cfg: LakeConfig,
) -> LakeGrid:
    """
    Create the full LakeGrid object.

    Parameters
    ----------
    cfg:
        LakeConfig instance.

    Returns
    -------
    grid:
        LakeGrid object containing all grid-related data.
    """
    x, y, X, Y, r, dx, dy, h = create_cartesian_grid(cfg)

    mask = create_disk_mask(r, cfg)

    boundary_band_topological = compute_topological_boundary_band(mask)

    (
        dist_to_boundary,
        boundary_band_1h,
        boundary_band_2h,
        boundary_band_3h,
    ) = compute_geometric_boundary_bands(
        r=r,
        mask=mask,
        h=h,
        cfg=cfg,
    )

    boundary_band = boundary_band_3h

    interior_unknowns = mask & (~boundary_band_topological)

    b, b_safe = compute_bathymetry(
        r=r,
        mask=mask,
        cfg=cfg,
    )

    psi0 = create_initial_stream_function(
        X=X,
        Y=Y,
        mask=mask,
        cfg=cfg,
    )

    return LakeGrid(
        radius=cfg.radius,
        x=x,
        y=y,
        X=X,
        Y=Y,
        r=r,
        dx=dx,
        dy=dy,
        h=h,
        mask=mask,
        dist_to_boundary=dist_to_boundary,
        boundary_band_topological=boundary_band_topological,
        boundary_band_1h=boundary_band_1h,
        boundary_band_2h=boundary_band_2h,
        boundary_band_3h=boundary_band_3h,
        boundary_band=boundary_band,
        interior_unknowns=interior_unknowns,
        b=b,
        b_safe=b_safe,
        psi0=psi0,
    )


def support_to_boundary_distance(
    field: np.ndarray,
    grid: LakeGrid,
    threshold: float = 1.0e-3,
) -> float:
    """
    Compute the distance from the significant support of a field to the boundary.

    Significant support is defined as:

        |field| > threshold * max_mask |field|.

    For a disk of radius R:

        dist(x, boundary) = R - r.

    Parameters
    ----------
    field:
        Scalar field.

    grid:
        LakeGrid instance.

    threshold:
        Relative threshold for significant support.

    Returns
    -------
    distance:
        Minimum distance to the boundary over the significant support.
        Returns NaN if the support is empty.
    """
    field_abs = np.abs(field)

    if not np.any(grid.mask):
        return np.nan

    max_value = np.max(field_abs[grid.mask])

    if max_value <= 0.0:
        return np.nan

    significant = grid.mask & (field_abs > threshold * max_value)

    if not np.any(significant):
        return np.nan

    return float(np.min(grid.radius - grid.r[significant]))


def print_grid_summary(
    grid: LakeGrid,
    cfg: LakeConfig,
) -> None:
    """
    Print a compact summary of the grid and bathymetry.

    Parameters
    ----------
    grid:
        LakeGrid instance.

    cfg:
        LakeConfig instance.
    """
    print("=" * 70)
    print("Grid and bathymetry summary")
    print("=" * 70)

    print(f"N                         = {cfg.N}")
    print(f"radius                    = {grid.radius:.8e}")
    print(f"dx                        = {grid.dx:.8e}")
    print(f"dy                        = {grid.dy:.8e}")
    print(f"h                         = {grid.h:.8e}")
    print(f"active grid points         = {int(grid.mask.sum())}")
    print(f"topological boundary nodes = {int(grid.boundary_band_topological.sum())}")
    print(f"boundary band 1h nodes     = {int(grid.boundary_band_1h.sum())}")
    print(f"boundary band 2h nodes     = {int(grid.boundary_band_2h.sum())}")
    print(f"boundary band 3h nodes     = {int(grid.boundary_band_3h.sum())}")

    print()
    print("Bathymetry")
    print(f"alpha                     = {cfg.alpha}")
    print(f"min b inside mask          = {np.min(grid.b[grid.mask]):.8e}")
    print(f"max b inside mask          = {np.max(grid.b[grid.mask]):.8e}")

    print()
    print("Initial stream function")
    dist_psi = support_to_boundary_distance(
        grid.psi0,
        grid,
        threshold=1.0e-6,
    )
    print(f"support distance psi0      = {dist_psi:.8e}")

    print("=" * 70)