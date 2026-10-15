# Reproduce figures from Section 4.3
include(joinpath(@__DIR__, "..", "Pi_Plot_Functions.jl"))
include(joinpath(@__DIR__, "..", "RK_Plot_Functions.jl"))
include(joinpath(@__DIR__, "..", "Pi_Methods_Data.jl"))

# ============================================================
# Setup Parameters
# ============================================================

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

# Cluster spectrum with r_fast = 5
S = make_eig_clusters(
    α = α,
    ε = ε,
    kappa_slow = 0.005,
    kappa_fast = 0.005
)
plotdir = joinpath(@__DIR__, "..", "Plots", "Pi")

# ============================================================
# Figure 4.4: Stability regions of PFEK2 and PRK4K2
# ============================================================

# Compute the maximal stable outer step sizes by bisection
dt_max_PFE, _ = find_dt_max_cluster_PI(
    K = 2,
    dt_inner = dt_inner,
    S = S,
    c_out = cout_PFE,
    A_out = Aout_PFE,
    b_out = bout_PFE,
    dt_max_search = 5.0,
    stab_tol = 0
);

dt_max_PRK4, _ = find_dt_max_cluster_PI(
    K = 2,
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
        K = 2,
        dt_inner = dt_inner,
        dt = dt_max_PFE,
    ),
    S = S,
    outdir = plotdir,
    filename = "PFE_cluster_r=5.png",
    rads = (slow = 1.3, fast = 30.0),
    stride = 1,
    spectrum_mode = :unscaled
)

plot_PI_slow_fast(
    cfg = (
        name = "PRK4",
        c = cout_RK4,
        A = Aout_RK4,
        b = bout_RK4,
        K = 2,
        dt_inner = dt_inner,
        dt = dt_max_PRK4,
    ),
    S = S,
    outdir = plotdir,
    filename = "PRK4_cluster_r=5.png",
    rads = (slow = 1.3, fast = 30.0),
    stride = 1,
    spectrum_mode = :unscaled
)

# ============================================================
# Figure 4.5: Effect of K on the PFE stability region
# ============================================================

plot_PI_slow_fast(
    cfg = (
        name = "PFE K=1",
        c = cout_PFE,
        A = Aout_PFE,
        b = bout_PFE,
        K = 1,
        dt_inner = dt_inner,
        dt = dt_max_PFE,
    ),
    S = S,
    outdir = plotdir,
    which = :fast,
    rads = (slow = nothing, fast = 25.0),
    stride = 1,
    filename = "PFE_cluster_r=5_K1.png",
    spectrum_mode = :unscaled
)

plot_PI_slow_fast(
    cfg = (
        name = "PFE K=1",
        c = cout_PFE,
        A = Aout_PFE,
        b = bout_PFE,
        K = 1,
        dt_inner = dt_inner,
        dt = dt_max_PFE,
    ),
    S = S,
    outdir = plotdir,
    which = :fast,
    rads = (slow = nothing, fast = 5.0),
    stride = 1,
    filename = "PFE_cluster_r=5_K1_zoomed.png",
    spectrum_mode = :unscaled
)


plot_PI_slow_fast(
    cfg = (
        name = "PFE K=2",
        c = cout_PFE,
        A = Aout_PFE,
        b = bout_PFE,
        K = 2,
        dt_inner = dt_inner,
        dt = dt_max_PFE,
    ),
    S = S,
    outdir = plotdir,
    which = :fast,
    rads = (slow = nothing, fast = 25.0),
    stride = 1,
    filename = "PFE_cluster_r=5_K2.png",
    spectrum_mode = :unscaled
)

# ============================================================
# Figure 4.6: Work-precision diagram for PFEK2 and PRK4K2
# ============================================================

p = work_precision_PI(
    problem = :cluster,
    S = S,
    T = 202,
    K = 2,
    dt_inner = dt_inner,
    annotate = true,
    outdir = plotdir,
    filename = "work_precision_PI_r=5_abs.png"
)