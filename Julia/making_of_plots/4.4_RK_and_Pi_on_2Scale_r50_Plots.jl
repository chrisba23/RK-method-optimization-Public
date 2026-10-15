# Reproduce figures from Section 4.4
include(joinpath(@__DIR__, "..", "Pi_Plot_Functions.jl"))
include(joinpath(@__DIR__, "..", "RK_Plot_Functions.jl"))
include(joinpath(@__DIR__, "..", "Pi_Methods_Data.jl"))

# 4-stage, fourth-order RK method as outer method
Aout_RK4 = [0.0  0.0  0.0  0.0;
            1/3  0.0  0.0  0.0;
           -1/3  1.0  0.0  0.0;
            1.0 -1.0  1.0  0.0]

bout_RK4 = [1/8, 3/8, 3/8, 1/8]
cout_RK4 = [0.0, 1/3, 2/3, 1.0]

# Forward Euler method as outer method
Aout_PFE = [0.0;;]
bout_PFE = [1.0]
cout_PFE = [0.0]

α = 1.0
ε = 1e-3
dt_inner = ε
K = 3

# Cluster spectrum with r_fast = 50
S = make_eig_clusters(
    α = α,
    ε = ε,
    kappa_slow = 0.05,
    kappa_fast = 0.05
)

pi_plotdir = joinpath(@__DIR__, "..", "Plots", "Pi")
rk_plotdir = joinpath(@__DIR__, "..", "Plots", "stability_regions")
# ============================================================
# Figure 4.7: Stability regions of PFEK3 and PRK4K3
# ============================================================

# Compute the maximal stable outer step sizes by bisection
dt_max_PFE, _ = find_dt_max_cluster_PI(
    K = K,
    dt_inner = dt_inner,
    S = S,
    c_out = cout_PFE,
    A_out = Aout_PFE,
    b_out = bout_PFE,
    dt_max_search = 5.0,
    stab_tol = 0
);

dt_max_PRK4, _ = find_dt_max_cluster_PI(
    K = K,
    dt_inner = dt_inner,
    S = S,
    c_out = cout_RK4,
    A_out = Aout_RK4,
    b_out = bout_RK4,
    dt_max_search = 5.0,
    stab_tol = 0
);


plot_PI_slow_fast(
    cfg = (
        name = "PFE",
        c = cout_PFE,
        A = Aout_PFE,
        b = bout_PFE,
        K = K,
        dt_inner = dt_inner,
        dt = dt_max_PFE,
    ),
    S = S,
    outdir = pi_plotdir,
    filename = "PFE_cluster_r=50_K=3.png",
    rads = (slow = 2.0, fast = 140.0),
    stride = 1,
    spectrum_mode = :unscaled
)

plot_PI_slow_fast(
    cfg = (
        name = "PRK4",
        c = cout_RK4,
        A = Aout_RK4,
        b = bout_RK4,
        K = K,
        dt_inner = dt_inner,
        dt = dt_max_PRK4,
    ),
    S = S,
    outdir = pi_plotdir,
    filename = "PRK4_cluster_r=50_K=3.png",
    rads = (slow = 2.0, fast = 140.0),
    stride = 1,
    spectrum_mode = :unscaled
)


# ============================================================
# Figure 4.8: Rejected 3Sstar p=2,s=8 method with Chebyshev basis
# ============================================================

matdir = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebychev", "2-Scale,r=50")
method_file = "3Sstar-2-8_2026-05-12T13-18-55.txt"


plot_compare_stability_methods(
    matdir,
    method_file;
    rad = nothing,
    rad_factor = 2.9,
    outdir = rk_plotdir,
    filename = "3Sstar_chebyshev_p2_s8_r50_compare.png"
)

# ============================================================
# Figure 4.9: Rejected 3Sstar p=3,s=10 method with Chebyshev basis
# ============================================================

method_file = "3Sstar-3-10_2026-04-22T22-48-46.txt"
plot_pseudocolor_comparison(
    matdir,
    method_file;
    rad = nothing,
    rad_factor = 1.8,
    outdir = rk_plotdir,
    filename = "3Sstar_chebyshev_p3_s10_r50_pseudocolor.png"
)





# ============================================================
# Figure 4.10: Rejected 2N p=2,s=9 method with Chebyshev basis
# ============================================================

matdir = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_2N", "Test_2N_chebychev", "2-Scale,r=50")
method_file = "2N-2-9_2026-04-28T06-14-44.txt"

plot_compare_stability_methods(
    matdir,
    method_file;
    rad = nothing,
    rad_factor = 2.9,
    outdir = rk_plotdir,
    filename = "2N_chebyshev_p2_s9_r50_compare.png"
)





# ============================================================
# Figure 4.11: 3Sstar p=4,s=14 method with Chebyshev interpolation
# using s - p - 1 matching conditions
# ============================================================
matdir = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=50","3Sstar_fast_focused")
method_file = "3Sstar-4-14_2026-04-18T15-07-52.txt"

plot_compare_stability_methods(
    matdir,
    method_file;
    rad = nothing,
    rad_factor = 1.3,
    outdir = rk_plotdir,
    filename = "3Sstar_chebyinterp_p4_s14_r50_compare.png"
)

# ============================================================
# Figure 4.12: Work-precision diagram for PFEK3 and PRK4K3
# ============================================================

p = work_precision_PI(
    problem = :cluster,
    S = S,
    T = 100,
    K = 3,
    dt_inner = dt_inner,
    annotate = true,
    outdir = pi_plotdir,
    filename = "work_precision_PI_r=50_abs.png"
)

