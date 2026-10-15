# Optimization of Explicit Runge-Kutta Methods for Multiscale Problems

This repository contains the code for the Master's thesis
*Optimization of Explicit Runge-Kutta Methods for Multiscale Problems*
by Christian Bartash at JGU Mainz.

## What this code does

The code constructs stability-optimized explicit Runge-Kutta methods for
multiscale spectra and compares them with the Projective Integration methods
PFE and PRK4 on the test problems of the thesis.

The workflow has two steps:

1. **MATLAB** — optimize the stability polynomial with PolyOpt and construct a
   corresponding low-storage Runge-Kutta method with RK-Opt. Each method is
   stored in a `.txt` file containing its Butcher and low-storage coefficients.
2. **Julia** — read the generated `.txt` files, apply the methods to the test
   problems, and compute the stability limits, error tables, and figures from
   Chapter 4. The Projective Integration methods are implemented directly in
   Julia.

## Requirements

**MATLAB**
* MATLAB Optimization Toolbox — required for `fmincon`.
* MATLAB Global Optimization Toolbox — required by RK-Opt.
* [CVX](https://cvxr.com/cvx) with the solvers SeDuMi and SDPT3 — required by PolyOpt.
  After installing CVX, run `cvx_setup` once.
* PolyOpt and RK-Opt are part of the
  [RK-Opt package](https://github.com/ketch/RK-Opt) by D. Ketcheson. 
  For the thesis, two modified versions of this package were used locally 
  but are not included in this public repository:
  * `RK-Opt-add_2N_fixed` — includes interpolation-based matching conditions,
    a grid-based line search over the step size in PolyOpt, and a low-storage
    stability evaluation in RK-Opt.
  * `RK-Opt-add_2N_fixed-derivative_conditions` — same as above, except that the order
    conditions for the Chebyshev interpolation basis are imposed directly through
    derivatives of the interpolation polynomial (Sec. 3.1.6).

  The remaining functionality is unchanged from the original package. Each
  MATLAB script adds the required version to the MATLAB path automatically,
  so no separate installation is needed.

**Julia 1.11.5**
Start Julia 1.11.5 from the root directory of the repository. Then run the
following commands once after cloning the repository:
```julia
import Pkg
Pkg.activate("Julia")
Pkg.instantiate()
```

The Julia code also uses the private `2025_hyp_RK_dev` repository by H. Ranocha to read in the optimized Runge-Kutta methods. 
This repository is not included in the public version.

## Optimization setups

Each spectrum considered in the thesis has its own directory in `Matlab/`.
The corresponding `Test_*.m` script specifies the order `p`, the number of
stages `s`, and the step-size search range. It calls PolyOpt and RK-Opt and
stores the resulting method in the same directory.

| Thesis section | Spectrum | Folder |
|---|---|---|
| 4.1 | isolated, λ = −1, −1000 | `Matlab/two_scale_lamS=-1_lamF=-1000/.../2-Scale,r=0` |
| 4.3 | cluster, r_fast = 5 | `Matlab/two_scale_lamS=-1_lamF=-1000/.../2-Scale,r=5` |
| 4.4 | cluster, r_fast = 50 | `Matlab/two_scale_lamS=-1_lamF=-1000/.../2-Scale,r=50` |
| 4.5 | moment model, ε = 1e-3 | `Matlab/hyp_upwind_model/lamS=-299.928_lamF=-1299.93_R=299.928` |
| 4.5 | moment model, ε = 1e-4 | `Matlab/hyp_upwind_model/lamS=-299.928_lamF=-10299.9_R=299.928` |
| 4.6 | generalized two-scale system | `Matlab/two_scale_lamS=-20_lamF=-2000` |

The subdirectories specify the low-storage class and the polynomial basis:
`Test_3Sstar_chebychev`, `Test_3Sstar_chebyinterpol`,
`Test_2N_chebychev`, and `Test_2N_chebychev_monotone`.

For the Chebyshev interpolation basis, the directory name also specifies the
distribution of the interpolation conditions: `fast_focused`,
`full_interval`, or `cluster_split`.

RK-Opt uses randomly generated initial guesses. Therefore, repeated runs with
the same setup may return different methods or may return no method.

## Reproducing the results of the thesis

**MATLAB** produces the PolyOpt tables of Sections 4.1 and 4.3:

| Script | Tables |
|---|---|
| `Matlab/making_of_data/4.1_r=0/Run_M_Dependence.m` | 4.2–4.5 |
| `Matlab/making_of_data/4.1_r=0/Run_linesearch.m` | 4.6 |
| `Matlab/making_of_data/4.3_r=5/Run_M_Dependence.m` | 4.13–4.16 |
| `Matlab/making_of_data/4.3_r=5/Run_linesearch.m` | 4.17 |

Here `4.1_r=0` refers to the isolated spectrum of Section 4.1 and `4.3_r=5` to the
cluster spectrum of Section 4.3. Results are written to `results_M_dependence/` and
`results_linesearch/` as `.xlsx`.

**Julia** produces all remaining tables and all figures. There is one script per section:

For example, the results for the isolated two-scale spectrum can be generated
from the root directory of the repository by running:

```julia
include("Julia/making_of_data/4.1_RK_and_Pi_on_2Scale_r0_Tables.jl")   # tables
include("Julia/making_of_plots/4.1_RK_and_Pi_on_2Scale_r0_Plots.jl")   # figures
```

and analogously for `4.3`, `4.4`, `4.5`, `4.6`. These scripts compute the maximal stable
outer step sizes of PFE/PRK4 for each K, read the optimized Runge-Kutta methods and check
their stability, and compute the error and work-precision results. Figures are written to
`Julia/Plots/`.

## References

* [RK-Opt](https://github.com/ketch/RK-Opt) — D. Ketcheson et al. The
  modified versions in `Extern-Repos/` are based on this package.
* [2025_hyp_RK_dev](https://github.com/ranocha/2025_hyp_RK_dev) — H. Ranocha
  (*private repository*).

See the thesis for the full bibliography of the methods used here.
