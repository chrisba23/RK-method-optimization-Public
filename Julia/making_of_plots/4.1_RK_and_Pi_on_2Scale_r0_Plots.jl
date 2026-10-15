# Reproduce figures from Section 4.1
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
K = 1

# Isolated two-scale spectrum Λ₀ = {-1, -1000}
S = make_eig_clusters(
    α = α,
    ε = ε,
    kappa_slow = 0.0,
    kappa_fast = 0.0
)

plotdir = joinpath(@__DIR__, "..", "Plots")
pi_plotdir = joinpath(plotdir, "Pi")

# ============================================================
# Figure 4.1: Exact solution of the two-scale problem
# ============================================================
plotdir = joinpath(@__DIR__, "..", "Plots")
p_full = plot_exact_two_scale(
    T = 5.0,
    outdir = plotdir,
    filename = "two_scale_exact_solution.png",
    legend_pos = :topright,
)

p_zoom = plot_exact_two_scale(
    T = 0.01,
    outdir = plotdir,
    filename = "two_scale_exact_solution_zoomed.png",
    legend_pos = :bottomright,
)



# ============================================================
# Figure 4.2: Stability regions of PFE and PRK4
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
)

dt_max_PRK4, _ = find_dt_max_cluster_PI(
    K = K,
    dt_inner = dt_inner,
    S = S,
    c_out = cout_RK4,
    A_out = Aout_RK4,
    b_out = bout_RK4,
    dt_max_search = 5.0,
    stab_tol = 0
)


# Save the legend used for the PI stability plots
save_PI_legend(
    outdir = pi_plotdir,
    filename = "PI_legend.pdf",
    fontsize = 10,
    ncols = 3,
    size = (800, 55)
)


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
    filename = "PFE_two_point_r=0_K=1.png",
    rads = (slow = 1.3, fast = 1.0),
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
    filename = "PRK4_two_point_r=0_K=1.png",
    rads = (slow = 1.3, fast = 1.0),
    spectrum_mode = :unscaled
)


# ============================================================
# Figure 4.3: Work-precision diagram for PFE and PRK4
# ============================================================

p = work_precision_PI(
    S = S,
    T = 100.0,
    K = 1,
    dt_inner = dt_inner,
    annotate = true,
    outdir = pi_plotdir,
    filename = "work_precision_PI_r=0_abs.png"
)

