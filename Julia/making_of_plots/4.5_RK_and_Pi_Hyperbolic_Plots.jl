# Reproduce figures from Section 4.5
include(joinpath(@__DIR__, "..", "Pi_Plot_Functions.jl"))
include(joinpath(@__DIR__, "..", "RK_Plot_Functions.jl"))
include(joinpath(@__DIR__, "..", "Pi_Methods_Data.jl"))

# ============================================================
# Setup Parameters
# ============================================================

matrixdir = joinpath(@__DIR__, "..", "Test_Problem_Matrices")

A_eps1e3 = readdlm(
    joinpath(matrixdir, "matrix_eps10e-3.csv"),
    ',',
    Float64,
)

A_eps1e4 = readdlm(
    joinpath(matrixdir, "matrix_eps10e-4.csv"),
    ',',
    Float64,
)

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


S1 = ComplexF64.(eigvals(A_eps1e3))
S2 = ComplexF64.(eigvals(A_eps1e4))

pi_plotdir = joinpath(@__DIR__, "..", "Plots", "Pi")
rk_plotdir = joinpath(@__DIR__, "..", "Plots", "stability_regions")

# ============================================================
# Figure 4.13: Stability regions for the semi-discrete moment
# model with ε = 1e-3 of PFEK1  and PRK4K2 
# ============================================================

λfast_eps1e3 = -1299.928
dt_inner_eps1e3 = -1 / λfast_eps1e3

dt_max_PFE_eps1e3, _ = find_dt_max_cluster_PI(
    K = 1,
    dt_inner = dt_inner_eps1e3,
    S = S1,
    c_out = cout_PFE,
    A_out = Aout_PFE,
    b_out = bout_PFE,
    dt_max_search = 0.05,
    stab_tol = 1e-12,
);

dt_max_PRK4_eps1e3, _ = find_dt_max_cluster_PI(
    K = 2,
    dt_inner = dt_inner_eps1e3,
    S = S1,
    c_out = cout_RK4,
    A_out = Aout_RK4,
    b_out = bout_RK4,
    dt_max_search = 0.05,
    stab_tol = 1e-12
);

plot_PI_full_spectrum(
    cfg = (
        name = "PFE",
        c = cout_PFE,
        A = Aout_PFE,
        b = bout_PFE,
        K = 1,
        dt_inner = dt_inner_eps1e3,
        dt = dt_max_PFE_eps1e3,
    ),
    S = S1,
    outdir = pi_plotdir,
    filename = "PFE_hyperbolic_K1_full.png",
    rad = 1300.0,
    stride = 1,
    spectrum_mode = :unscaled
)

plot_PI_full_spectrum(
    cfg = (
        name = "PRK4",
        c = cout_RK4,
        A = Aout_RK4,
        b = bout_RK4,
        K = 2,
        dt_inner = dt_inner_eps1e3,
        dt = dt_max_PRK4_eps1e3,
    ),
    S = S1,
    outdir = pi_plotdir,
    filename = "PRK4_hyperbolic_K2_full.png",
    rad = 1300.0,
    stride = 1,
    spectrum_mode = :unscaled
)

# ============================================================
# Figure 4.14: Stability regions for the semi-discrete moment
# model with ε = 1e-4 of PFEK1  and PRK4K1
# ============================================================

λfast_eps1e4 = -10299.928
dt_inner_eps1e4 = -1 / λfast_eps1e4


dt_max_PFE_eps1e4, _ = find_dt_max_cluster_PI(
    K = 1,
    dt_inner = dt_inner_eps1e4,
    S = S2,
    c_out = cout_PFE,
    A_out = Aout_PFE,
    b_out = bout_PFE,
    dt_max_search = 0.02,
    stab_tol = 1e-12,
);

dt_max_PRK4_eps1e4, _ = find_dt_max_cluster_PI(
    K = 1,
    dt_inner = dt_inner_eps1e4,
    S = S2,
    c_out = cout_RK4,
    A_out = Aout_RK4,
    b_out = bout_RK4,
    dt_max_search = 0.02,
    stab_tol = 1e-12,
);

plot_PI_slow_fast(
    cfg = (
        name = "PFE",
        c = cout_PFE,
        A = Aout_PFE,
        b = bout_PFE,
        K = 1,
        dt_inner = dt_inner_eps1e4,
        dt = dt_max_PFE_eps1e4,
    ),
    S = S2,
    outdir = pi_plotdir,
    filename = "PFE_hyperbolic_K1.png",
    rads = (slow = 800.0, fast = 900.0),
    stride = 1,
    spectrum_mode = :unscaled
)

plot_PI_slow_fast(
    cfg = (
        name = "PRK4",
        c = cout_RK4,
        A = Aout_RK4,
        b = bout_RK4,
        K = 1,
        dt_inner = dt_inner_eps1e4,
        dt = dt_max_PRK4_eps1e4,
    ),
    S = S2,
    outdir = pi_plotdir,
    filename = "PRK4_hyperbolic_K1.png",
    rads = (slow = 800.0, fast = 900.0),
    stride = 1,
    spectrum_mode = :unscaled
)




# ============================================================
# Figure 4.15: Comparison of stability regions for the semi-discrete moment
# for chebyshev Basis with p=4 and s=16 
# ============================================================
matdir = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyshev")
method_file = "3Sstar-4-16_2026-05-12T22-49-21.txt"

plot_compare_stability_methods(
    matdir,
    method_file;
    rad = nothing,
    rad_factor = 2.9,
    outdir = rk_plotdir,
    filename = "3Sstar_chebyshev_p4_s16_hyperbolic_compare.png"
)

