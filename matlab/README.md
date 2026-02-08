# MATLAB solver for spherical-cap snap-through (axisymmetric BVP)

Run in MATLAB:

```matlab
results = solve_spherical_cap_snapthrough();
```

## What is implemented
- The solver uses the **same non-dimensional equations in `\phi`** as in your manuscript (including full `tan(\phi)`, `cot(\phi)`, `sin(\phi)`, `cos(\phi)` terms; no shallow-angle replacement inside the governing equations).
- Unknown state is `[\bar w, \bar w', \bar w'', \bar N_\varphi, \bar N_\theta, q]`, where `q` is the transverse bracket term used to keep a first-order system for `bvp4c`.
- Continuation between load steps uses `bvpinit(sol.x, @(phi) deval(sol,phi))` so the guess is always vector-valued.

## Load model options
Set at the top of the script:
- `p.loadModel = 'point'` (default): concentrated apex load enforced through apex BC in `q` form.
- `p.loadModel = 'patch'`: optional regularized pressure patch over `0 <= phi <= phiPatch`.
- `p.patchNormalization` controls how `Fbar` is converted to patch pressure in `loadModel='patch'`:
  - `'domain'` (default): equivalent load is matched on the solved interval `[epsPole, phiPatch]` (best comparison with point mode in this code).
  - `'fullcap'`: equivalent load is matched on `[0, phiPatch]`; the omitted `0..epsPole` core is added back through a pole `q` jump correction.
- `p.patchLoadSign`: sign of pressure forcing in the transverse `q'` equation. Use `+1` (default) to match the concentrated-load convention used in this script; if your historical sign convention is opposite, set `-1`.

## Control modes
- `p.controlMode = 'load'` (default): prescribe `FbarList`, solve for deflection.
- `p.controlMode = 'displacement'`: prescribe `p.wApexList`, solve for force `Fbar` such that `wbar(eps)\approx wbar(0)` matches each target.

## User parameters
- `p.nu`: Poisson ratio.
- `p.gamma`: slenderness ratio `R/h`.
- `p.phi0`: cap half-angle in radians.
- `p.epsPole`: small positive start angle to avoid coordinate singularity at exactly `phi=0`.
- `p.phiPatch`: patch half-angle (used only for `loadModel='patch'`).
- `FbarList`: load continuation values, `Fbar = F*R/D` (used in load control).
- `p.wApexList`: prescribed apex displacement targets for displacement control.
- `p.forcePerturb`: finite-difference force increment used to estimate `dw/dF` during displacement control Newton updates.
- `p.dispTol`: displacement-control matching tolerance on `wbar(eps)-wTarget`.
- `p.maxDispNewton`: maximum Newton iterations per prescribed displacement target.
- `p.maxLoadStep`: automatic sub-stepping size in load continuation (smaller = more robust at high load).
- `p.residualNodes`: number of sample points used to print a posteriori ODE/BC residual indicators after each step.
- `p.residualTrim`: fraction trimmed near both boundaries when evaluating ODE residual indicator (reduces pole singular-coordinate contamination).
- `p.residualQuantile`: robust quantile (default 0.95) used for ODE residual reporting.

## Output
- `results(i).Fbar`: load level `\mathcal F = FR/D`.
- `results(i).wApex`: approximate apex deflection `\bar w(0)` sampled at `phi=epsPole`.
- `results(i).sol`: full `bvp4c` solution structure for post-processing.

## Why point and patch can differ
- Even with equal total `Fbar`, a true point load and a finite patch are different physical load distributions, so responses are only expected to coincide asymptotically as `phiPatch -> 0`.
- To reduce discrepancy with minimal edits, keep `loadModel='patch'`, set `patchNormalization='domain'`, and use small `phiPatch` (e.g. `0.005`–`0.02` depending on mesh/conditioning).

## Meaning of the second plotted curve
- The second figure is **not** a second load-deflection branch.
- It is the spatial deflection profile `\bar w(\phi)` at the **last load step** (`results(end).Fbar`).
- X-axis is position along the cap (`phi` from pole to edge), Y-axis is local non-dimensional displacement.
- This tells you *shape* at one load, while the first figure tells you *global response* across loads.

## High-load robustness (minimal strategy implemented)
- The solver now automatically splits large load jumps into sub-steps (`maxLoadStep`) to mimic continuation-relaxation and reduce solver stalls at high load.
- After each sub-step, it prints:
  - `ODE-relres-q95`: robust interior relative ODE residual indicator (95th percentile by default).
  - `ODE-relres-max`: maximum interior relative ODE residual (spike-sensitive, diagnostic only).
  - `BC-res-inf`: max boundary-condition residual.
- If you see warnings or long runtime, reduce `p.maxLoadStep` (e.g. from `10` to `2` or `1`).

## Displacement-control algorithm
- For each target `wTarget` in `p.wApexList`, the code iteratively updates `Fbar` so that `wbar(eps)-wTarget -> 0`.
- It uses a finite-difference Newton step: solve at `Fbar` and at `Fbar + dF`, estimate `dw/dF`, then update `Fbar <- Fbar - r/(dw/dF)`.
- The converged `Fbar` is the force that satisfies the apex displacement request.
