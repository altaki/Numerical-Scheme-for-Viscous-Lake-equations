# Degenerate Viscous Lake Equations — Numerical Discretization

This document describes **the continuous model**, **the choices of discretization** (space/time), and **the weighted projection** used to enforce the constraint \( \nabla\cdot(b\,u)=0 \) within the context of **viscous lake equations** with degenerate bathymetry \( b\ge0 \) that vanishes at the boundary.

The main references are:

- Al Taki (2017), *Viscosity effect on the degenerate lake equations*, **Nonlinear Analysis** 148, 30–60 — existence/uniqueness in weighted spaces and power-type weights. [doi:10.1016/j.na.2016.09.017][AT2017]
- Al Taki & Lacave (2023), *Degenerate lake equations: classical solutions and vanishing viscosity limit*, **Nonlinearity** 36, 653–678 — classical solutions of the degenerate inviscid model and **limit \( \nu\to0 \)** in the Navier case. [doi:10.1088/1361-6544/aca865][AT-L2023]
- Jiu, Niu, Wu (2012), *Vanishing viscosity limits for the degenerate lake equations with Navier boundary conditions*, **Nonlinearity** 25, 641–655 — Navier framework and inviscid limit. [PDF IOP][JNW2012]

https://doi.org/10.1016/j.na.2016.09.017  
[AT-L2023]: https://doi.org/10.1088/1361-6544/aca865  
https://iopscience.iop.org/article/10.1088/0951-7715/25/3/641/pdf

> **Note** — The numerical setup below follows the spirit of weak forms and weighted constraints used in the literature (conservative discretization and weighted projection). For more analytical and physical details on Navier, see also Iftimie–Sueur (2010) and Wang–Wang–Xin (2010).

---

## 1) Continuous Model (recap)

Let \( \Omega\subset\mathbb{R}^2 \) be bounded, and a bathymetry \( b(x)\ge0 \) (e.g., **Muckenhoupt-type weight**) that vanishes at the boundary: typically
$$
 b(x) = \operatorname{dist}(x,\partial\Omega)^{\alpha},\qquad \alpha>0.
$$
The weighted viscous system (common 2D form) writes
\[
\partial_t( b\,u ) + \nabla\cdot( b\,u\otimes u )
- 2\mu\,\nabla\cdot\Big( b\big(D(u)+ (\nabla\cdot u)I\big) \Big)
+ b\,\nabla p = 0,\qquad \nabla\cdot(b\,u)=0,
\]
where \( D(u)=\tfrac12(\nabla u+\nabla u^T) \). In the limit \( \mu\to0 \) we recover the inviscid model (“lake/Euler anelastic”), studied in [AT‑L2023] and [JNW2012].

### Boundary Conditions (Navier)
On \( \partial\Omega \):
- **Impermeability**: \( (b\,u)\cdot n = 0 \).
- **Navier slip**: \( 2b\,(D(u)\,n)\cdot\tau + \eta\,b\,(u\cdot\tau)=0 \) (\( \eta\ge0 \)).

These conditions avoid Dirichlet-type boundary layers and are standard for justifying the **limit \( \mu\to0 \)** (see [JNW2012], [AT‑L2023]).

---

## 2) Discretization (Cartesian grid, conservative form + weighted projection)

### 2.1 Grid and Indexing
- Domain \( [0,L_x]\times[0,L_y] \), mesh \( n_x\times n_y \).
- **Cell-centered variables**: \( u=(u_x,u_y),\ b,\ p \) tabulated on indices \( (i,j) \) with \( i=0..n_x-1 \) for **x** and \( j=0..n_y-1 \) for **y** (i.e. `indexing='ij'`).
- **Face weights** (arithmetic mean):  
  \( b_{i+\frac12,j}=\tfrac12(b_{i+1,j}+b_{i,j}) \),  
  \( b_{i,j+\frac12}=\tfrac12(b_{i,j+1}+b_{i,j}) \).

### 2.2 Time Scheme (Fractional Step)
Let \( \Delta t \) be the time step. One step \( n\to n+1 \):

1. **Explicit Advection (conservative weighted form)**. For each component \( u_k \):
   \[
     u_k^{\*} = u_k^n - \Delta t\,\frac{1}{b}\,\nabla\cdot\big( b\,u^n\, u_k^n\big)
                 + \Delta t\,\text{Diff}(u_k^n).
   \]
   **Fluxes** at faces are **upwinded**:  
   \( F^{x}_{i+\frac12,j}=b_{i+\frac12,j}\,(u_x)_{i+\frac12,j}\,(u_k)^{\text{up}}_{i+\frac12,j} \),  
   \( F^{y}_{i,j+\frac12}=b_{i,j+\frac12}\,(u_y)_{i,j+\frac12}\,(u_k)^{\text{up}}_{i,j+\frac12} \),  
   and \( \nabla\cdot(\cdot) \) is **flux divergence** (conservative finite differences).

2. **Weighted Projection** (enforce \( \nabla\cdot(b u^{n+1})=0 \)). Solve
   \[
     \nabla\cdot(b\,\nabla\phi) = \frac{1}{\Delta t}\,\nabla\cdot(b\,u^{\*}),
     \qquad \text{Neumann BC: } \partial_n\phi=0.
   \]
   Then **correct without \( b \) factor**:  
   \( u^{n+1} = u^{\*} - \nabla\phi \) and set \( p^{n+1}=\phi \) (up to a factor).

3. **Boundary Conditions**. Default:
   - **Normal**: \( u\cdot n=0 \) (if \( b>0 \) at the boundary; otherwise trivially satisfied if \( b=0 \)).
   - **Tangential**: **free-Navier** version \( (\eta=0) \) via \( \partial_n(u\cdot\tau)=0 \).
     Robin \( (\eta>0) \) is implemented in ghost cells:  
     \( 2b(D(u)n)\cdot\tau + \eta b(u\cdot\tau)=0 \).

> **Why this projection?** It is the Euler-Lagrange of  
> \( \min_v \tfrac12\int b|v-u^{\*}|^2 \) s.t. \( \nabla\cdot(b v)=0 \),  
> hence \( u^{n+1}=u^{\*}-\nabla\phi \) and the weighted elliptic above (see [JNW2012]).

### 2.3 Diffusion (two variants)
- **Simple variant (per component)**: \( \text{Diff}(u)=\nu\,\frac{1}{b}\,\nabla\cdot(b\,\nabla u) \).
  5-point discretization **with variable coefficients**:
  \[
  (\nabla\cdot(b\,\nabla u))_{i,j} \approx \frac{b_{i+\frac12,j}(u_{i+1,j}-u_{i,j})-b_{i-\frac12,j}(u_{i,j}-u_{i-1,j})}{\Delta x^2}
   + \frac{b_{i,j+\frac12}(u_{i,j+1}-u_{i,j})-b_{i,j-\frac12}(u_{i,j}-u_{i,j-1})}{\Delta y^2}.
  \]
- **“Tensorial” variant**:  
  \( -2\mu\,\nabla\cdot( b(D(u)+(\nabla\cdot u)I)) \). Requires assembling a vector operator (mix of cross derivatives); preferable in **FEM**.

> **Important** — Avoid `np.roll` (implicit periodic conditions) if one wants to enforce Navier.

---

## 3) Stability & Practical Remarks
- **CFL** advection: \( \Delta t \lesssim C\,\min(\Delta x,\Delta y)/\|u\|_\infty \).
- **Degenerate weights**: the elliptic \( \nabla\cdot(b\nabla\phi) \) may be **ill-conditioned** near boundaries when \( b\to0 \) — prefer **harmonic averages** for faces and a **multigrid preconditioner** (AMG).
- **Weighted energy**: without forcing, \( \tfrac12\int b|u|^2 \) should decrease (diffusion) — good consistency test.

---

## 4) Minimal Validation
1. **Weighted divergence**: check \( \|\nabla\cdot(bu^{n+1})\|_{L^2} \) after projection.
2. **Smooth case with \( b\equiv\text{const} \)**: compare to 2D Navier–Stokes (reduces to classical case).
3. **Grid/time convergence**: decrease of residuals and error upon refinement.
4. **Study \( \nu\to0 \)**: measure \( \|u_\nu-u_0\|_{L^2(b)} \) (inviscid reference computed with adapted scheme) as in [AT‑L2023], [JNW2012].

---

## 5) References
- B. Al Taki, *Viscosity effect on the degenerate lake equations*, **Nonlinear Analysis** 148 (2017), 30–60. [doi10.1016/j.na.2016.09.017
- B. Al Taki & C. Lacave, *Degenerate lake equations: classical solutions and vanishing viscosity limit*, **Nonlinearity** 36 (2023), 653–678. doi:10.1088/1361-6544/aca865
- Q. Jiu, D. Niu, J. Wu, *Vanishing viscosity limits for the degenerate lake equations with Navier boundary conditions*, **Nonlinearity** 25 (2012), 641–655. [IOP PDF](https://iopscience.iop.org/article/10.1088/0951-7715/25/3/641/pdf)
- D. Iftimie, F. Sueur, *Viscous boundary layers for Navier–Stokes with Navier slip*, preprint (2010). [PDF](https://math.univ-lyon1.fr/%7Eiftimie/ARTICLES/geominv.pdf)
- X.-P. Wang, Y.-G. Wang, Z. Xin, *Boundary layers in incompressible Navier–Stokes with Navier boundary conditions*, **Commun. Math. Sci.** 8(4) (2010), 965–998. [PDF](https://intlpress.com/site/pub/files/_fulltext/journals/cms/2010/0008/0004/CMS-2010-0008-0004-a010.pdf)
