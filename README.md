# Numerical-Scheme-for-Viscous-Lake-equations

# Équations des lacs visqueux (dégénérées) — Discrétisation numérique

Ce document décrit **le modèle continu**, **les choix de discrétisation** (espace/temps), et **la
projection pondérée** utilisée pour imposer la contrainte \(\nabla\cdot(b\,u)=0\) dans le cadre
des **équations des lacs visqueux** avec bathymétrie dégénérée \(b\ge0\) qui s’annule au bord.

Les références de base sont :

- Al Taki (2017), *Viscosity effect on the degenerate lake equations*, **Nonlinear Analysis** 148, 30–60 — existence/unicité en espaces pondérés et poids de type puissance. [doi:10.1016/j.na.2016.09.017][AT2017]
- Al Taki & Lacave (2023), *Degenerate lake equations: classical solutions and vanishing viscosity limit*, **Nonlinearity** 36, 653–678 — solutions classiques du modèle inviscide dégénéré et **limite \(\nu\to0\)** sous Navier. [doi:10.1088/1361-6544/aca865][AT-L2023]
- Jiu, Niu, Wu (2012), *Vanishing viscosity limits for the degenerate lake equations with Navier boundary conditions*, **Nonlinearity** 25, 641–655 — cadre Navier et limite inviscide. [PDF IOP][JNW2012]

[AT2017]: https://doi.org/10.1016/j.na.2016.09.017
[AT-L2023]: https://doi.org/10.1088/1361-6544/aca865
[JNW2012]: https://iopscience.iop.org/article/10.1088/0951-7715/25/3/641/pdf

> **Remarque** — La mise en place numérique ci‑dessous suit l’esprit des formes faibles et
> contraintes pondérées utilisées dans la littérature (discrétisation conservatrice et projection
> pondérée). Pour des détails analytiques et physiques supplémentaires sur Navier, voir aussi
> Iftimie–Sueur (2010) et Wang–Wang–Xin (2010). 

---

## 1) Modèle continu (rappel)

Soit \(\Omega\subset\mathbb{R}^2\) borné, et une bathymétrie \(b(x)\ge0\) (ex. **poids de type Muckenhoupt**)
qui s’annule au bord : typiquement
\[
 b(x) = \operatorname{dist}(x,\partial\Omega)^{\alpha},\qquad \alpha>0.
\]
Le système visqueux pondéré (forme courante en 2D) s’écrit
\[
\partial_t( b\,u ) + \nabla\cdot( b\,u\otimes u )
- 2\mu\,\nabla\cdot\Big( b\big(D(u)+ (\nabla\cdot u)I\big) \Big)
+ b\,\nabla p = 0,\qquad \nabla\cdot(b\,u)=0,
\]
où \(D(u)=\tfrac12(\nabla u+\nabla u^T)\). En limite \(\mu\to0\) on retrouve le modèle inviscide
(“lake/Euler anelastic”), étudié dans [AT‑L2023] et [JNW2012].

### Conditions au bord (Navier)
Sur \(\partial\Omega\) :
- **Imperméabilité** : \((b\,u)\cdot n = 0\).
- **Glissement de Navier** : \(2b\,(D(u)\,n)\cdot\tau + \eta\,b\,(u\cdot\tau)=0\) (\(\eta\ge0\)).

Ces conditions évitent les couches limites de type Dirichlet et sont standard pour justifier la
**limite \(\mu\to0\)** (voir [JNW2012], [AT‑L2023]).

---

## 2) Discrétisation (grille cartésienne, forme conservative + projection pondérée)

### 2.1 Grille et indexation
- Domaine \([0,L_x]\times[0,L_y]\), maillage \(n_x\times n_y\).
- Variables **au centre de cellule** : \(u=(u_x,u_y),\ b,\ p\) tabulés sur indices \((i,j)\) avec
  \(i=0..n_x-1\) pour **x** et \(j=0..n_y-1\) pour **y** (i.e. `indexing='ij'`).
- Poids **aux faces** (moyenne arithmétique) : 
  \(b_{i+\frac12,j}=\tfrac12(b_{i+1,j}+b_{i,j})\),
  \(b_{i,j+\frac12}=\tfrac12(b_{i,j+1}+b_{i,j})\).

### 2.2 Schéma en temps (fractionnaire)
Soit \(\Delta t\) le pas de temps. Un pas \(n\to n+1\) :

1. **Advection explicite (forme conservative pondérée)**. Pour chaque composante \(u_k\) :
   \[
     u_k^{\*} = u_k^n - \Delta t\,\frac{1}{b}\,\nabla\cdot\big( b\,u^n\, u_k^n\big)
                 + \Delta t\,\text{Diff}(u_k^n).
   \]
   Les **flux** aux faces sont **upwindés** :
   \(
   F^{x}_{i+\frac12,j}=b_{i+\frac12,j}\,(u_x)_{i+\frac12,j}\,(u_k)^{\text{up}}_{i+\frac12,j}
   \),
   \(
   F^{y}_{i,j+\frac12}=b_{i,j+\frac12}\,(u_y)_{i,j+\frac12}\,(u_k)^{\text{up}}_{i,j+\frac12}
   \),
   et \(\nabla\cdot(\cdot)\) est **la divergence de flux** (différences finies conservatrices).

2. **Projection pondérée** (impose \(\nabla\cdot(b u^{n+1})=0\)). Résoudre
   \[
     \nabla\cdot(b\,\nabla\phi) = \frac{1}{\Delta t}\,\nabla\cdot(b\,u^{\*}),
     \qquad \text{C.L. Neumann : } \partial_n\phi=0.
   \]
   Puis **corriger sans facteur \(b\)** :
   \(u^{n+1} = u^{\*} - \nabla\phi\) et poser \(p^{n+1}=\phi\) (à un facteur près).

3. **Conditions au bord**. Par défaut :
   - **Normal** : \(u\cdot n=0\) (si \(b>0\) au bord ; sinon satisfait trivialement si \(b=0\)).
   - **Tangentiel** : version **free‑Navier** \((\eta=0)\) via \(\partial_n(u\cdot\tau)=0\).
     Un Robin \((\eta>0)\) s’implémente en cellule fantôme :
     \(2b(D(u)n)\cdot\tau + \eta b(u\cdot\tau)=0\).

> **Pourquoi cette projection ?** C’est l’EL de
> \(\min_v \tfrac12\int b|v-u^{\*}|^2\) s.c. \(\nabla\cdot(b v)=0\), 
> d’où \(u^{n+1}=u^{\*}-\nabla\phi\) et l’elliptique pondéré ci‑dessus (voir [JNW2012]).

### 2.3 Diffusion (deux variantes)
- **Variante simple (par composante)** : \(\text{Diff}(u)=\nu\,\frac{1}{b}\,\nabla\cdot(b\,\nabla u)\).
  Discrétisation 5‑points **à coefficients variables** :
  \[
  (\nabla\cdot(b\,\nabla u))_{i,j} \approx \frac{b_{i+\frac12,j}(u_{i+1,j}-u_{i,j})-b_{i-\frac12,j}(u_{i,j}-u_{i-1,j})}{\Delta x^2}
   + \frac{b_{i,j+\frac12}(u_{i,j+1}-u_{i,j})-b_{i,j-\frac12}(u_{i,j}-u_{i,j-1})}{\Delta y^2}.
  \]
- **Variante “tensorielle”** fidèle :
  \(-2\mu\,\nabla\cdot( b(D(u)+(\nabla\cdot u)I))\). Demande d’assembler un opérateur vectoriel
  (mélange de dérivées croisées) ; préférable en **FEM**.

> **Important** — Éviter `np.roll` (conditions périodiques implicites) si l’on veut imposer Navier.

---

## 3) Stabilité & remarques pratiques
- **CFL** advection : \(\Delta t \lesssim C\,\min(\Delta x,\Delta y)/\|u\|_\infty\).
- **Poids dégénéré** : l’elliptique \(\nabla\cdot(b\nabla\phi)\) peut être **mal conditionné** près des
  bords lorsque \(b\to0\) — privilégier **moyennes harmoniques** pour les faces et un
  **préconditionneur multigrille** (AMG).
- **Énergie pondérée** : sans forçage, \(\tfrac12\int b|u|^2\) doit décroître (diffusion) — bon test de
  consistance.

---

## 4) Validation minimale
1. **Divergence pondérée** : contrôler \(\|\nabla\cdot(bu^{n+1})\|_{L^2}\) après projection.
2. **Cas lisse avec \(b\equiv\text{const}\)** : comparer à Navier–Stokes 2D (réduit au cas classique).
3. **Convergence maillage/temps** : décroissance des résidus et de l’erreur au raffinement.
4. **Étude \(\nu\to0\)** : mesurer \(\|u_\nu-u_0\|_{L^2(b)}\) (référence inviscide calculée avec schéma
   adapté) comme dans [AT‑L2023], [JNW2012].

---

## 5) Références
- B. Al Taki, *Viscosity effect on the degenerate lake equations*, **Nonlinear Analysis** 148 (2017), 30–60. [doi:10.1016/j.na.2016.09.017](https://doi.org/10.1016/j.na.2016.09.017)
- B. Al Taki & C. Lacave, *Degenerate lake equations: classical solutions and vanishing viscosity limit*, **Nonlinearity** 36 (2023), 653–678. [doi:10.1088/1361-6544/aca865](https://doi.org/10.1088/1361-6544/aca865)
- Q. Jiu, D. Niu, J. Wu, *Vanishing viscosity limits for the degenerate lake equations with Navier boundary conditions*, **Nonlinearity** 25 (2012), 641–655. [IOP PDF](https://iopscience.iop.org/article/10.1088/0951-7715/25/3/641/pdf)
- D. Iftimie, F. Sueur, *Viscous boundary layers for Navier–Stokes with Navier slip*, preprint (2010). [PDF](https://math.univ-lyon1.fr/%7Eiftimie/ARTICLES/geominv.pdf)
- X.-P. Wang, Y.-G. Wang, Z. Xin, *Boundary layers in incompressible Navier–Stokes with Navier boundary conditions*, **Commun. Math. Sci.** 8(4) (2010), 965–998. [PDF](https://intlpress.com/site/pub/files/_fulltext/journals/cms/2010/0008/0004/CMS-2010-0008-0004-a010.pdf)
