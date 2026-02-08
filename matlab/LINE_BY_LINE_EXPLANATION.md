# Line-by-line scientific explanation of `solve_spherical_cap_snapthrough.m`

This document explains the MATLAB solver line by line (grouped in contiguous line ranges for readability), with the physical and numerical meaning of each part.

---

## 1) Function header and mathematical scope (lines 1–24)

- **Line 1** defines the main MATLAB function `solve_spherical_cap_snapthrough()` that returns a `results` structure containing the continuation solutions.
- **Lines 2–4** state that the solver is for an axisymmetric spherical-cap BVP using your non-dimensional formulation in the meridional angle `phi`, explicitly avoiding hidden shallow-angle substitution.
- **Lines 6–12** define the unknown state vector
  \[
  y = [\bar w,\, \bar w',\, \bar w'',\, \bar N_\varphi,\, \bar N_\theta,\, q]^T,
  \]
  where `q` is introduced as an auxiliary transverse resultant to reduce differential order for `bvp4c`.
- **Lines 13–18** list the governing equations used in the code:
  1. membrane equilibrium,
  2. transverse equilibrium written in terms of `q`,
  3. compatibility equation with exact trigonometric terms.
- **Lines 20–23** describe apex-load handling:
  - default: true concentrated load (point mode) through a pole boundary condition on `q`,
  - optional: patch regularization mode for robustness checks.

Scientific role: these lines explicitly state the strong form model and clarify the numerical variable transformation from a high-order shell equation to a first-order BVP system.

---

## 2) User parameters and load continuation setup (lines 25–47)

- **Lines 25–29** define material/geometric nondimensional parameters:
  - `nu`: Poisson ratio,
  - `gamma = R/h`: shell slenderness,
  - `phi0`: cap edge angle,
  - `epsPole`: small positive lower domain boundary to avoid evaluating singular geometric coefficients exactly at `phi=0`.
- **Lines 31–33** define load-model control:
  - `p.loadModel = 'point'` enforces concentrated apex load through BC,
  - `p.phiPatch` is active only if `loadModel='patch'`.
- **Lines 35–36** define continuation loads `FbarList`; each entry is a separate nonlinear BVP solve with previous solution used as initial guess.
- **Lines 38–41** build initial mesh and initial guess object for `bvp4c`:
  - line 39 sets initial `p.Fbar` for consistency in first guess generation,
  - line 41 uses callback form `@(phi) guess_profile(phi,p)`.
- **Lines 43–45** preallocate output structure fields (`Fbar`, `sol`, `wApex`) to avoid dynamic resizing.
- **Line 47** defines solver tolerances (`RelTol`, `AbsTol`) and maximum mesh refinement (`NMax`).

Scientific role: continuation is critical for nonlinear shell problems with possible snap-through branches; stepping load gradually improves convergence and branch tracking.

---

## 3) Nonlinear continuation loop (lines 49–65)

- **Line 49** starts loop over prescribed dimensionless loads.
- **Line 50** updates current load value in parameter struct.
- **Lines 52–54** solve the BVP at current load using:
  - ODE residual function `ode_system`,
  - boundary residual function `bc_system`,
  - previous/initial guess `solinit`.
- **Line 57** updates continuation seed with
  `bvpinit(sol.x, @(phi) deval(sol,phi))`,
  ensuring vector-valued guess function and avoiding matrix-guess compatibility errors in some MATLAB environments.
- **Line 59** evaluates solution near apex (`phi=epsPole`) to estimate apex displacement.
- **Lines 60–62** store current load and solution in `results`.
- **Line 64** prints convergence progress for each load step.

Scientific role: each converged step approximates a quasi-static equilibrium at a load level; the printed sequence helps identify possible sudden displacement changes associated with snap-through behavior.

---

## 4) Post-processing plots (lines 67–83)

- **Lines 67–73** create load-deflection curve:
  - x-axis: `Fbar = FR/D`,
  - y-axis: apex deflection estimate `wbar(eps)`.
- **Lines 75–82** plot full deflection profile \(\bar w(\phi)\) at the highest load.

Scientific role: load-deflection response is the primary indicator of stiffness loss, limit points, and post-buckling/snap-through branches.

---

## 5) Initial guess function `guess_profile` (lines 85–103)

- **Line 85** begins local function used by `bvpinit`.
- **Lines 87–90** construct a smooth small-deflection guess for \(\bar w\), \(\bar w'\), \(\bar w''\) using normalized coordinate `s` from apex to edge.
- **Lines 92–93** initialize membrane forces to zero.
- **Lines 95–99** initialize `q` depending on load model:
  - in point mode, start from physically consistent apex-load value,
  - in patch mode, start from zero.
- **Lines 101–102** assemble vector and force column shape with `y(:)`.

Scientific role: a physically consistent initial guess significantly stabilizes nonlinear BVP convergence, especially under high geometric nonlinearity.

---

## 6) ODE system function `ode_system` (lines 105–152)

### 6.1 State unpacking and geometric coefficients (lines 107–120)

- **Lines 107–115** unpack state variables and parameters.
- **Lines 117–120** compute exact trigonometric factors `sin`, `cos`, `tan`, `cot` needed by spherical geometry operators.

Scientific role: keeping exact trigonometric operators preserves full spherical geometry (non-shallow form).

### 6.2 Reconstruction of third derivative from auxiliary variable `q` (lines 122–124)

- `q` definition is algebraically rearranged to compute `wppp = d^3 w/dphi^3`.

Scientific role: this converts the original higher-order transverse equation into first-order system form acceptable by `bvp4c`.

### 6.3 Membrane equilibrium equation (lines 126–127)

- Implements
  \[
  (N_\varphi\sin\phi)' - N_\theta\cos\phi = 0
  \Rightarrow
  N_\varphi' = (N_\theta-N_\varphi)\cot\phi.
  \]

Scientific role: enforces in-surface force balance along meridian/circumference coupling.

### 6.4 Compatibility equation solved for \(N_\theta'\) (lines 129–133)

- `G` stores non-derivative terms from compatibility relation.
- Then
  \[
  N_\theta' = \nu N_\varphi' - \frac{G}{\tan\phi}.
  \]

Scientific role: closes membrane part by coupling strains and curvature-induced geometric nonlinearity via \(\bar w\), \(\bar w'\).

### 6.5 Transverse equilibrium with selectable load model (lines 135–149)

- **Point mode (line 138):**
  \[
  q' = \gamma(N_\varphi+N_\theta)\sin\phi,
  \]
  with concentrated-load effect imposed via pole BC, not distributed source.
- **Patch mode (lines 140–147):**
  computes equivalent constant pressure in a small cap region so integrated load equals `Fbar`, then subtracts distributed forcing term.
- **Lines 148–149** guard against unsupported modes.

Scientific role: separates mathematical model choices for singular loading treatment while keeping same shell equations.

### 6.6 Return first-order derivative vector (line 151)

- Returns
  \([w', w'', w''', N_\varphi', N_\theta', q']^T\).

---

## 7) Boundary-condition function `bc_system` (lines 154–176)

### 7.1 Apex target for `q` by load model (lines 159–166)

- **Point mode:**
  \[
  q(\epsilon) \approx -\gamma\mathcal{F}/(2\pi),
  \]
  implementing concentrated apex load in weak-limit form.
- **Patch mode:** `q` pole target is zero because load is distributed in-domain.

### 7.2 Full BC residual vector (lines 168–175)

Six BC equations for six first-order unknowns:
1. `ya(2)=0`: \(w'(\epsilon)=0\) (pole symmetry),
2. `ya(4)-ya(5)=0`: \(N_\varphi=N_\theta\) at pole (regularity),
3. `ya(6)-qPoleTarget=0`: apex-load condition in `q` form,
4. `yb(1)=0`: clamped edge deflection,
5. `yb(2)=0`: clamped edge slope,
6. `yb(5)-nu*yb(4)=0`: edge membrane relation from your boundary setup.

Scientific role: BC set encodes both geometric constraints (clamp/symmetry) and load application mechanism.

---

## 8) What this solver does **not** do by itself

- It does not automatically branch-switch at limit points; it follows the branch reachable by continuation in `FbarList`.
- It does not include additional constitutive terms for explicit Gurtin–Murdoch surface elasticity yet; current implementation corresponds to the bulk-shell nondimensional equations you provided.
- It does not evaluate exactly at `phi=0`; it uses `epsPole>0` and asymptotic BC consistency.

---

## 9) Minimal usage recipe

```matlab
results = solve_spherical_cap_snapthrough();
```

Then inspect:

```matlab
F = [results.Fbar];
wA = [results.wApex];
plot(F, wA, 'o-'); grid on;
```

If you want strict concentrated-load treatment, keep `p.loadModel = 'point'`.
Use `'patch'` only for regularized sensitivity checks.
