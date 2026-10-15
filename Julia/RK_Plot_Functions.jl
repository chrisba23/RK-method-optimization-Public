# ============================================================
# Runge-Kutta stability and solution plots
#
# Loads PolyOpt plot data, evaluates the stability
# polynomial, compares different stability-function evaluations,
# and creates plots for the clusters and ODE solutions.
# ============================================================

include(joinpath(@__DIR__, "setup.jl"))

# Plot settings
const STABLE_COLOR   = RGB(0.15, 0.75, 0.75)
const SPECTRUM_COLOR = :red
const BOUNDARY_COLOR = :black
stable_cgrad() = cgrad([STABLE_COLOR, STABLE_COLOR])


function save_plot(plt, outdir, method_file; filename = nothing)
    mkpath(outdir)

    if filename === nothing
        base, _ = splitext(method_file)
        plotfile = joinpath(outdir, base * ".png")
    else
        plotfile = joinpath(outdir, filename)
    end

    savefig(plt, plotfile)
    println("Saved plot: ", abspath(plotfile))

    return plotfile
end

# Build a plotting grid around the scaled fast cluster.
function fast_cluster_grid(S, h_star;
    rad = nothing,
    rad_factor = 1.5,
    n = 600,
    stride = 12)

    _, S_fast = split_slow_fast(S)

    z = h_star .* S_fast

    z0 = mean(z)
    r = maximum(abs.(z .- z0))
    rad_use = rad === nothing ? rad_factor * max(r, 1.0) : rad

    xs = range(real(z0) - rad_use, real(z0) + rad_use, length = n)
    ys = range(imag(z0) - rad_use, imag(z0) + rad_use, length = n)
    Z = [x + im * y for y in ys, x in xs]

    z_spec = z[1:stride:end]

    return z0, rad_use, xs, ys, Z, z_spec
end

# Evaluate the RK stability function in power-series, Butcher,
# or low-storage form.
function stability_values(mode, alg, method_type, A, b, Z)
    if mode == "powers"
        coeffs = stability_poly_coeffs(A, b)
        R = map(z -> evalpoly(z, coeffs), Z)
        title = "Power-series"

    elseif mode == "butcher"
        R = R_grid(A, b, Z)
        title = "Butcher"

    elseif mode == "lowstorage"
        R = R_lowstorage(alg, method_type, Z)
        title = "Low-storage"

    else
        error("Unknown mode: $mode")
    end

    return R, title
end

# Plot the RK stability region and PolyOpt boundary around the
# slow and fast clusters.
function plot_slow_fast_clusters(matdir, method_file;
    rads = (slow = nothing, fast = nothing),
    rad_factor = 1.5,
    n = 600,
    stride = 12,
    outdir = joinpath(@__DIR__, "plots", "stability_regions"),
    filename = nothing)

    experiment = load_experiment(matdir, method_file)

    method_type = String(experiment["method_type"])
    alg = load_method(matdir, method_file, method_type)

    A = alg.coefficients.A
    b = alg.coefficients.b

    h_star = Float64(experiment["h_star"])
    S = vec(experiment["S"])

    S_slow, S_fast = split_slow_fast(S)

    z_slow = h_star .* S_slow
    z_fast = h_star .* S_fast

    clusters = [z_slow, z_fast]
    titles = ["Slow cluster", "Fast cluster"]

    centers = [mean(z_slow), mean(z_fast)]

    auto_rads = [
        rad_factor * max(maximum(abs.(z_slow .- centers[1])), 1.0),
        rad_factor * max(maximum(abs.(z_fast .- centers[2])), 1.0)
    ]

    plot_rads = [
        rads.slow === nothing ? auto_rads[1] : rads.slow,
        rads.fast === nothing ? auto_rads[2] : rads.fast
    ]

    plots = Plots.Plot[]

    for j in [2, 1]
        z_spec = clusters[j]
        z0 = centers[j]
        rad = plot_rads[j]

        z_plot = z_spec[1:stride:end]

        xs = range(real(z0) - rad, real(z0) + rad, length = n)
        ys = range(imag(z0) - rad, imag(z0) + rad, length = n)

        Z = [x + im*y for y in ys, x in xs]

        R_rk = abs.(R_grid(A, b, Z))
        R_poly = abs.(eval_polyopt_reference(experiment, Z))

        stable_fill = ifelse.(R_rk .<= 1.0, 1.0, NaN)

        plt = Plots.heatmap(
            xs, ys, stable_fill;
            color = stable_cgrad(),
            colorbar = false,
            aspect_ratio = :equal,
            xlabel = "Re(z)",
            ylabel = "Im(z)",
            titlefontsize = 16,
            guidefontsize = 14,
            tickfontsize = 12,
            title = titles[j],
            legend = false,
            grid = true,
            framestyle = :box
        )

        Plots.contour!(
            plt,
            xs, ys, R_poly;
            levels = [1.0],
            color = BOUNDARY_COLOR,
            linewidth = 2,
            label = false
        )

        Plots.scatter!(
            plt,
            real.(z_plot),
            imag.(z_plot);
            color = SPECTRUM_COLOR,
            markersize = 3,
            markerstrokewidth = 0,
            label = false
        )

        Plots.xlims!(plt, real(z0) - rad, real(z0) + rad)
        Plots.ylims!(plt, imag(z0) - rad, imag(z0) + rad)

        push!(plots, plt)
    end

    leg = Plots.plot(
        framestyle = :none,
        grid = false,
        xlims = (0, 1),
        ylims = (0, 1),
        xticks = false,
        yticks = false,
        widen = false,
        legend = :top,
        legend_columns = 3,
        legendfontsize = 9
    )

    Plots.plot!(
        leg,
        [NaN], [NaN];
        color = STABLE_COLOR,
        linewidth = 8,
        label = "RK stable region"
    )

    Plots.plot!(
        leg,
        [NaN], [NaN];
        color = BOUNDARY_COLOR,
        linewidth = 2,
        label = "PolyOpt boundary"
    )

    Plots.scatter!(
        leg,
        [NaN], [NaN];
        color = SPECTRUM_COLOR,
        markersize = 3,
        markerstrokewidth = 0,
        label = "Scaled spectrum"
    )

    layout = Plots.@layout [
        p1 p2
        l{0.10h}
    ]

    plt_all = Plots.plot(
        plots[1],
        plots[2],
        leg;
        layout = layout,
        size = (1050, 540),
        margin = 3Plots.mm,
        bottom_margin = 8Plots.mm,
        left_margin = 7Plots.mm
    )

    save_plot(
        plt_all,
        outdir,
        method_file;
        filename = filename
    )
    return plt_all
end


# Plot the exact solution of the two-scale ODE test problem.
function plot_exact_two_scale(;
    T = 5.0,
    α = 1.0,
    ε = 1e-3,
    u10 = 1.0,
    u20 = 0.0,
    n_points = 2000,
    legend_pos = :topright,
    outdir = joinpath(@__DIR__, "plots"),
    filename = "exact_two_scale.png"
)

    tvals = range(0.0, T; length = n_points)

    vals = [
        exact_two_scale(
            t;
            α = α,
            ε = ε,
            u10 = u10,
            u20 = u20,
        )
        for t in tvals
    ]

    u1 = first.(vals)
    u2 = last.(vals)

    plt = Plots.plot(
        tvals,
        u2;
        label = "u₂(t)",
        xlabel = "t",
        ylabel = "uᵢ(t)",
        linewidth = 4.5,
        legend = legend_pos,
        size = (900, 550),
        dpi = 300,
        guidefontsize = 20,
        tickfontsize = 15,
        legendfontsize = 16,
        left_margin = 7Plots.mm,
        bottom_margin = 6Plots.mm,
        right_margin = 3Plots.mm,
        top_margin = 3Plots.mm,
    )

    Plots.plot!(
        plt,
        tvals,
        u1;
        label = "u₁(t)",
        linewidth = 4.5,
        linestyle = :dash,
    )

    mkpath(outdir)

    plotfile = joinpath(outdir, filename)

    savefig(plt, plotfile)
    println("Saved plot: ", abspath(plotfile))

    return plt
end


# Compare the RK stability regions obtained from the low-storage,
# Butcher, and power-series evaluations near the fast cluster.
function plot_compare_stability_methods(matdir, method_file;
    rad = nothing,
    rad_factor = 1.5,
    n = 600,
    outdir = joinpath(@__DIR__, "plots", "stability_regions"),
    filename = nothing)

    experiment = load_experiment(matdir, method_file)
    method_type = String(experiment["method_type"])
    alg = load_method(matdir, method_file, method_type)

    A = alg.coefficients.A
    b = alg.coefficients.b

    h_star = Float64(experiment["h_star"])
    S = vec(experiment["S"])

    z0, rad, xs, ys, Z, z_spec = fast_cluster_grid(S, h_star; rad = rad, rad_factor = rad_factor, n = n)

    R_polyopt = abs.(eval_polyopt_reference(experiment, Z))

    modes = ["lowstorage", "butcher", "powers"]

    R_methods = []
    titles = String[]

    for mode in modes
        R, title = stability_values(mode, alg, method_type, A, b, Z)
        push!(R_methods, abs.(R))
        push!(titles, title)
    end

    plots = Plots.Plot[]

    for j in 1:3
        stable_fill = ifelse.(R_methods[j] .<= 1.0, 1.0, NaN)

        plt = Plots.heatmap(
            xs, ys, stable_fill;
            color = stable_cgrad(),
            colorbar = false,
            aspect_ratio = :equal,
            xlabel = "Re(z)",
            ylabel = "Im(z)",
            guidefontsize = 11,
            tickfontsize = 9,
            title = titles[j],
            legend = false,
            grid = true,
            framestyle = :box
        )

        Plots.contour!(
            plt,
            xs, ys, R_polyopt;
            levels = [1.0],
            color = BOUNDARY_COLOR,
            linewidth = 2,
            label = false
        )

        Plots.scatter!(
            plt,
            real.(z_spec),
            imag.(z_spec);
            color = SPECTRUM_COLOR,
            markersize = 3,
            markerstrokewidth = 0,
            label = false
        )
        Plots.xlims!(plt, real(z0) - rad, real(z0) + rad)
        Plots.ylims!(plt, imag(z0) - rad, imag(z0) + rad)

        push!(plots, plt)
    end

    leg = Plots.plot(
        framestyle = :none,
        grid = false,
        xlims = (0, 1),
        ylims = (0, 1),
        xticks = false,
        yticks = false,
        widen = false,
        legend = :top,
        legend_columns = 3,
        legendfontsize = 11
    )

    Plots.plot!(
        leg, [NaN], [NaN];
        color = STABLE_COLOR,
        linewidth = 8,
        label = "RK stable region"
    )

    Plots.plot!(
        leg, [NaN], [NaN];
        color = BOUNDARY_COLOR,
        linewidth = 2,
        label = "PolyOpt boundary"
    )

    Plots.scatter!(
        leg, [NaN], [NaN];
        color = SPECTRUM_COLOR,
        markersize = 3,
        markerstrokewidth = 0,
        label = "Scaled spectrum"
    )

    layout = Plots.@layout [
        p1 p2 p3
        l{0.12h}
    ]

    plt_all = Plots.plot(
        plots[1],
        plots[2],
        plots[3],
        leg;
        layout = layout,
        size = (1350, 500),
        margin = 1Plots.mm,
        legend_columns = 3
    )

    save_plot(
        plt_all,
        outdir,
        method_file;
        filename = filename
    )
    return plt_all
end

# Compare the stability-function values obtained from the low-storage,
# Butcher, and power-series evaluations using pseudocolor plots.
function plot_pseudocolor_comparison(matdir, method_file;
    rad = nothing,
    rad_factor = 1.5,
    n = 600,
    outdir = joinpath(@__DIR__, "plots", "stability_regions"),
    filename = nothing)

    experiment = load_experiment(matdir, method_file)
    method_type = String(experiment["method_type"])
    alg = load_method(matdir, method_file, method_type)

    A = alg.coefficients.A
    b = alg.coefficients.b
    h_star = Float64(experiment["h_star"])
    S = vec(experiment["S"])

    z0, rad, xs, ys, Z, z_spec = fast_cluster_grid(
        S,
        h_star;
        rad = rad,
        rad_factor = rad_factor,
        n = n
    )

    modes = ["lowstorage", "butcher", "powers"]

    levels = [0.0, 0.25, 0.5, 0.75, 0.9, 1.0, 1.05, 1.25]

    colors = [
        RGB(0.15, 0.20, 0.65),
        RGB(0.25, 0.50, 0.85),
        RGB(0.30, 0.75, 0.70),
        RGB(0.75, 0.85, 0.45),
        RGB(0.95, 0.90, 0.35),
        RGB(0.85, 0.85, 0.85),
        RGB(0.65, 0.65, 0.65)
    ]

    clipmax = 1.25
    cgrad_R = cgrad(
        colors,
        length(colors),
        categorical = true
    )

    plots = Plots.Plot[]

    for (j, mode) in enumerate(modes)
        R_sel, plot_title = stability_values(
            mode,
            alg,
            method_type,
            A,
            b,
            Z
        )

        V = abs.(R_sel)
        Vplot = min.(V, clipmax)

        plt = Plots.contourf(
            xs,
            ys,
            Vplot;
            levels = levels,
            color = cgrad_R,
            colorbar = false,
            aspect_ratio = :equal,
            xlabel = "Re(z)",
            ylabel = j == 1 ? "Im(z)" : "",
            title = plot_title,
            guidefontsize = 11,
            tickfontsize = 9,
            clims = (0, clipmax),
            grid = true,
            framestyle = :box,
            legend = false
        )

        Plots.contour!(
            plt,
            xs,
            ys,
            V;
            levels = [1.0],
            color = BOUNDARY_COLOR,
            linewidth = 2,
            label = false
        )

        Plots.scatter!(
            plt,
            real.(z_spec),
            imag.(z_spec);
            color = SPECTRUM_COLOR,
            markerstrokecolor = :black,
            markerstrokewidth = 0.5,
            markersize = 3,
            label = false
        )

        Plots.xlims!(
            plt,
            real(z0) - rad,
            real(z0) + rad
        )

        Plots.ylims!(
            plt,
            imag(z0) - rad,
            imag(z0) + rad
        )

        push!(plots, plt)
    end
    cbar = Plots.plot(
        framestyle = :none,
        grid = false,
        xlims = (0, 1.25),
        ylims = (0, 1.25),
        xticks = false,
        yticks = false,
        widen = false,
        legend = false,
        title = "|R_RK(z)|",
        titlefontsize = 9,
        margin = 0Plots.mm
    )

    x0, x1 = 0.18, 0.48

    for k in eachindex(colors)
        Plots.plot!(
            cbar,
            Plots.Shape(
                [x0, x1, x1, x0],
                [
                    levels[k],
                    levels[k],
                    levels[k + 1],
                    levels[k + 1]
                ]
            );
            color = colors[k],
            linecolor = colors[k],
            label = false
        )
    end

    Plots.plot!(
        cbar,
        [x0, x1],
        [1.0, 1.0];
        color = BOUNDARY_COLOR,
        linewidth = 2,
        label = false
    )

    for lev in [0.0, 0.25, 0.5, 0.75, 1.0, 1.25]
        lab = lev == 0.0 ? "0" : string(round(lev, digits = 2))

        Plots.annotate!(
            cbar,
            x1 + 0.10,
            lev,
            Plots.text(lab, 8, :left)
        )
    end

    # Legend
    leg = Plots.plot(
        framestyle = :none,
        grid = false,
        xlims = (0, 1),
        ylims = (0, 1),
        xticks = false,
        yticks = false,
        widen = false,
        legend = :top,
        legend_columns = 2,
        legendfontsize = 9
    )

    Plots.plot!(
        leg,
        [NaN],
        [NaN];
        color = BOUNDARY_COLOR,
        linewidth = 2,
        label = "|R_RK(z)| = 1"
    )

    Plots.scatter!(
        leg,
        [NaN],
        [NaN];
        color = SPECTRUM_COLOR,
        markerstrokecolor = :black,
        markerstrokewidth = 0.5,
        markersize = 4,
        label = "Scaled spectrum"
    )

    layout = Plots.@layout [
        p1 p2 p3 c{0.06w}
        l{0.09h}
    ]

    plt_all = Plots.plot(
        plots[1],
        plots[2],
        plots[3],
        cbar,
        leg;
        layout = layout,
        size = (1420, 500),
        margin = 1Plots.mm,
        bottom_margin = 8Plots.mm,
        top_margin = 4Plots.mm,
        left_margin = 5Plots.mm,
        right_margin = 5Plots.mm,
        legend_columns = 2
    )

    save_plot(
        plt_all,
        outdir,
        method_file;
        filename = filename
    )

    return plt_all
end
