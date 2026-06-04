# Numerical Scheme for Viscous Lake Equations

This repository contains a **masked Cartesian prototype** for the numerical study of degenerate lake equations and a qualitative investigation of the vanishing-viscosity limit.

$$
\Omega = \lbrace{ (x,y)\in\mathbb{R}^2 : x^2+y^2<1 \rbrace}.
$$

It is approximated numerically by a masked Cartesian grid on the square

$$
[-1,1]^2.
$$

The bathymetry is

$$
b(x,y) = (1-r^2)^\alpha,
\qquad
r=\sqrt{x^2+y^2},
\qquad
\alpha=0.4.
$$

Near the boundary,

$$
1-r^2=(1-r)(1+r)\sim 2\,\mathrm{dist}(x,\partial\Omega),
$$

so

$$
b(x,y)\sim \mathrm{dist}(x,\partial\Omega)^\alpha.
$$

Thus \(b\to0\) at the shore.

The default configuration uses

$$
0<\alpha<\frac 12,
$$

which is the degenerate regime considered in the motivating mathematical framework.

---

## 2. Inviscid lake equations

The inviscid problem is solved in vorticity formulation:

$$
\partial_t\omega + u\cdot\nabla\omega = 0.
$$

The velocity is reconstructed from a stream function \(\psi\) through the singular elliptic problem

$$
\mathrm{div}\left(\frac1b\nabla\psi\right)=b\omega.
$$

The velocity is then recovered by

$$
u=\frac1b\nabla^\perp\psi.
$$

The convention used throughout the code is

$$
\nabla^\perp\psi=(-\partial_y\psi,\partial_x\psi),
$$

and

$$
\mathrm{curl}u=\partial_xu_y-\partial_yu_x.
$$

Therefore,

$$
\omega=\frac{\mathrm{curl}u}{b}.
$$

The sparse elliptic reconstruction solves the SPD system

$$
A_{\rm sing}\psi=-b\omega,
$$

where

$$
A_{\rm sing}=-\mathrm{div}\left(\frac1b\nabla\right).
$$

The singular elliptic operator is implemented in `lake/elliptic.py`.

---

## 3. Viscous lake equations

The viscous problem is treated in velocity formulation.

The continuous model motivating the prototype is

$$
\partial_t(bu_\mu)
+\mathrm{div}(bu_\mu\otimes u_\mu)
-2\mu\mathrm{div}\left(bD(u_\mu)+b\mathrm{div}(u_\mu)I\right)
+b\nabla p_\mu=0,
$$

with

$$
\mathrm{div}(bu_\mu)=0.
$$

The symmetric gradient is

$$
D(u)=\frac12(\nabla u+\nabla u^T).
$$

The continuous theory includes Navier-type boundary conditions:

$$
bu_\mu\cdot n=0,
$$

and

$$
2b\left(D(u_\mu)\cdot n+\mathrm{div}(u_\mu)n\right)\cdot\tau
+\eta_\mu b(u_\mu\cdot\tau)=0.
$$

However, this prototype does **not** impose the full Navier slip-with-friction boundary condition exactly. The boundary is represented by a Cartesian mask, and boundary effects are only approximated and monitored through diagnostics.

The viscous solver is implemented in `lake/viscous.py` as a masked Cartesian predictor-projection prototype.

---

## 4. Vanishing-viscosity comparison

The vanishing-viscosity comparison is performed only at the velocity level.

The numerical error is measured in the weighted velocity norm

$$
\lVert u_\mu-u\rVert _{L^2_b}
= \left(\int_\Omega b |u_\mu-u|^2\,dx\right)^{1/2}.
$$

On the masked grid this is approximated by

$$
\lVert u_\mu-u\rVert _{L^2_b}=
\left(
\sum_{\Omega_h} b_{ij}|u_{\mu,ij}-u_{ij}|^2\,dx\,dy
\right)^{1/2}.
$$

The prototype does **not** compare viscous and inviscid vorticities.

A formal theoretical reference rate has the form

$$
O\left(\mu^{(1-\beta)/2}\right),
$$

under assumptions including a Navier friction scaling

$$
0\leq\eta_\mu\leq \eta\mu^{-\beta},
\qquad
\beta<1.
$$

In this prototype, the observed numerical error contains several contributions:

- viscosity error;
- spatial discretization error;
- time discretization error;
- elliptic solver error;
- projection consistency error;
- boundary-mask geometry error;
- projection regularization error.

Therefore, the fitted slope in \(\mu\) must be interpreted qualitatively unless all these error sources are controlled.

---

## 5. Project structure

```text
lake_vanishing_viscosity/
│
├── main_notebook.ipynb
│
├── lake/
│   ├── __init__.py
│   ├── config.py
│   ├── grid.py
│   ├── operators.py
│   ├── elliptic.py
│   ├── inviscid.py
│   ├── projection.py
│   ├── viscous.py
│   ├── diagnostics.py
│   └── plotting.py
│
├── scripts/
│   └── run_quick_test.py
│
├── README.md
└── .gitignore

The project was designed to separate clearly:

1. the **inviscid lake equations** in vorticity-stream formulation;
2. the **viscous lake equations** in velocity formulation;
3. the **vanishing-viscosity comparison** in the weighted velocity norm \(L^2_b\).

This is a prototype implementation. It is not a boundary-fitted FEM/FVM solver and should not be interpreted as a rigorous numerical verification of the theoretical vanishing-viscosity convergence rate.

---



