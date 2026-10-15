# Reproduce figures from Section 4.6
include(joinpath(@__DIR__, "..", "Pi_Plot_Functions.jl"))
include(joinpath(@__DIR__, "..", "RK_Plot_Functions.jl"))
include(joinpath(@__DIR__, "..", "Pi_Methods_Data.jl"))

# ============================================================
# Setup Parameters
# ============================================================

matrixdir = joinpath(@__DIR__, "..", "Test_Problem_Matrices")

A = readdlm(
    joinpath(matrixdir, "matrixN20x2.csv"),
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


S = ComplexF64.(eigvals(A))

λfast_center = -2000.0
dt_inner = -1 / λfast_center
K = 8
plotdir = joinpath(@__DIR__, "..", "Plots", "Pi")

# ============================================================
# Figure 4.16: Stability regions of PFEK8 and PRK4K8 for the
# generalized two-scale system
# ============================================================

dt_max_PFE, _ = find_dt_max_cluster_PI(
    K = K,
    dt_inner = dt_inner,
    S = S,
    c_out = cout_PFE,
    A_out = Aout_PFE,
    b_out = bout_PFE,
    dt_max_search = 0.2,
    stab_tol = 0,
);

dt_max_PRK4, _ = find_dt_max_cluster_PI(
    K = K,
    dt_inner = dt_inner,
    S = S,
    c_out = cout_RK4,
    A_out = Aout_RK4,
    b_out = bout_RK4,
    dt_max_search = 0.2,
    stab_tol = 0,
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
    outdir = plotdir,
    filename = "PFE_generalisation_K=8.png",
    rads = (slow = 30.0, fast = 1340.0),
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
    outdir = plotdir,
    filename = "PRK4_generalisation_K=8.png",
    rads = (slow = 30.0, fast = 1340.0),
    stride = 1,
    spectrum_mode = :unscaled
)

# ============================================================
# Pseudospectrum and stability-region plots 
# ============================================================

function plot_pseudospectrum_slow_fast(;
    A_system,
    S,
    cfg,
    method_type,
    outdir,
    filename,
    eta_ps = 40.0,
    rads = (slow = 75.0, fast = 1300.0),
    n = 150
)
    if method_type == :PI
        λ_inner = cfg.dt_inner / cfg.dt
        _, A_method, b_method = create_PRK_tableau(
            cfg.c, cfg.A, cfg.b, cfg.K, λ_inner
        )
    else
        A_method = cfg.A
        b_method = cfg.b
    end

    S_slow, S_fast = split_slow_fast(S)
    clusters = [S_slow, S_fast]
    centers = mean.(clusters)
    plot_rads = [rads.slow, rads.fast]

    plots = Plots.Plot[]

    for j in [2, 1]
        spec = clusters[j]
        center = centers[j]
        rad = plot_rads[j]

        xs = range(real(center) - rad, real(center) + rad, length = n)
        ys = range(imag(center) - rad, imag(center) + rad, length = n)
        Z = [x + im*y for y in ys, x in xs]

        N_system = size(A_system, 1)
        I_system = Matrix{ComplexF64}(I, N_system, N_system)
        A_complex = ComplexF64.(A_system)

        resolvent_norm = reshape(
            [
                opnorm((z * I_system - A_complex) \ I_system,2)
                for z in vec(Z)
            ],
            size(Z)
        )

        R = abs.(R_grid(A_method, b_method, cfg.dt .* Z))
        stable_fill = ifelse.(R .<= 1.0, 1.0, NaN)

        plt = Plots.heatmap(
            xs, ys, stable_fill;
            color = stable_cgrad(),
            colorbar = false,
            xlabel = "Re(λ)",
            ylabel = "Im(λ)",
            title = cfg.name *
                    (j == 1 ? ": slow region" : ": fast region"),
            aspect_ratio = :equal,
            guidefontsize = 11,
            tickfontsize = 9,
            legend = false,
            grid = true,
            framestyle = :box
        )

        Plots.contour!(
            plt, xs, ys, resolvent_norm;
            levels = [1 / eta_ps],
            color = :blue,
            linewidth = 2,
            colorbar = false,
            label = false
        )

        Plots.contour!(
            plt, xs, ys, R;
            levels = [1.0],
            color = BOUNDARY_COLOR,
            linewidth = 2,
            colorbar = false,
            label = false
        )

        Plots.scatter!(
            plt, real.(spec), imag.(spec);
            color = SPECTRUM_COLOR,
            markersize = 3,
            markerstrokewidth = 0,
            label = false
        )
        Plots.xlims!(plt, real(center) - rad, real(center) + rad)
        Plots.ylims!(plt, imag(center) - rad, imag(center) + rad)
        push!(plots, plt)
    end

    plt = Plots.plot(
        plots...;
        layout = (1, 2),
        size = (1100, 440),
        margin = 2Plots.mm,
        bottom_margin = 8Plots.mm,
        left_margin = 5Plots.mm,
        right_margin = 5Plots.mm
    )

    mkpath(outdir)
    plotfile = joinpath(outdir, filename)
    savefig(plt, plotfile)
    println("Saved plot: ", abspath(plotfile))

    return plt
end

# ============================================================
# Figure 4.17: Pseudospectrum and stability-region plots 

# PFEK8
plot_pseudospectrum_slow_fast(
    A_system = A,
    S = S,
    cfg = (
        name = "PFEK8",
        c = cout_PFE,
        A = Aout_PFE,
        b = bout_PFE,
        K = K,
        dt_inner = dt_inner,
        dt = dt_max_PFE
    ),
    method_type = :PI,
    outdir = plotdir,
    filename = "PFE_generalisation_K8_pseudospectrum.png",
    eta_ps = 40.0,
    rads = (slow = 75.0, fast = 1300.0),
    n = 100
)


# PRK4K8
plot_pseudospectrum_slow_fast(
    A_system = A,
    S = S,
    cfg = (
        name = "PRK4K8",
        c = cout_RK4,
        A = Aout_RK4,
        b = bout_RK4,
        K = K,
        dt_inner = dt_inner,
        dt = dt_max_PRK4
    ),
    method_type = :PI,
    outdir = plotdir,
    filename = "PRK4_generalisation_K8_pseudospectrum.png",
    eta_ps = 40.0,
    rads = (slow = 75.0, fast = 1300.0),
    n = 100
)


# Optimized p=3,s=18 Runge-Kutta method using Chebyshev basis
matdir_3S_chebyshev = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_3Sstar", "Test_3Sstar_chebychev")
method_file_RK = "3Sstar-3-18_2026-06-28T04-25-30.txt"
experiment_RK = load_experiment(
    matdir_3S_chebyshev,
    method_file_RK
)
alg_RK = load_method(
    matdir_3S_chebyshev,
    method_file_RK,
    "3Sstar"
)
dt_RK = Float64(experiment_RK["h_star"])


plot_pseudospectrum_slow_fast(
    A_system = A,
    S = S,
    cfg = (
        name = "RK p=3, s=18",
        A = alg_RK.coefficients.A,
        b = alg_RK.coefficients.b,
        dt = dt_RK
    ),
    method_type = :RK,
    outdir = plotdir,
    filename = "RK_p3_s18_generalisation_pseudospectrum.png",
    eta_ps = 40.0,
    rads = (slow = 75.0, fast = 1300.0),
    n = 100
)

# ============================================================
# save legend for pseudospectrum plots
# ============================================================

function save_pseudospectrum_legend(;
    outdir,
    filename = "pseudospectrum_legend.pdf",
    fontsize = 10,
    size = (900, 55)
)
    leg = Plots.plot(
        framestyle = :none,
        grid = false,
        xlims = (0, 1),
        ylims = (0, 1),
        xticks = false,
        yticks = false,
        widen = false,
        legend = :top,
        legend_columns = 4,
        legendfontsize = fontsize,
        size = size,
        margin = 0Plots.mm,
        left_margin = 0Plots.mm,
        right_margin = 0Plots.mm,
        top_margin = 0Plots.mm,
        bottom_margin = 0Plots.mm
    )

    Plots.plot!(
        leg, [NaN], [NaN];
        color = STABLE_COLOR,
        linewidth = 7,
        label = "Stable region"
    )

    Plots.plot!(
        leg, [NaN], [NaN];
        color = BOUNDARY_COLOR,
        linewidth = 2,
        label = "|R(Δt λ)| = 1"
    )

    Plots.plot!(
        leg, [NaN], [NaN];
        color = :blue,
        linewidth = 2,
        label = "η-pseudospectrum"
    )

    Plots.scatter!(
        leg, [NaN], [NaN];
        color = SPECTRUM_COLOR,
        markersize = 4,
        markerstrokewidth = 0,
        label = "Spectrum"
    )

    mkpath(outdir)
    legendfile = joinpath(outdir, filename)
    savefig(leg, legendfile)

    println("Saved legend: ", abspath(legendfile))

    return leg
end

save_pseudospectrum_legend(
    outdir = plotdir,
    filename = "pseudospectrum_legend.pdf"
)


# ============================================================
# Figure 4.18: Work-precision diagram for PFEK8 and PRK4K8
# ============================================================

u0 = vcat(collect(1.0:20), zeros(20))
p = work_precision_PI(
    problem = :generalized,
    S = S,
    A_sys = A,
    u0 = u0,
    T = 50.0,
    K = 8,
    dt_inner = dt_inner,
    factors = [1, 0.9, 0.8,0.6],
    annotate = true,
    outdir = plotdir,
    filename = "work_precision_generalized_K8.png",
)