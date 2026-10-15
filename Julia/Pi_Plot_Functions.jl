# ============================================================
# Projective integration stability-region plots
# ============================================================

include(joinpath(@__DIR__, "setup.jl"))

# shared Plot settings
const STABLE_COLOR   = RGB(0.15, 0.75, 0.75)
const SPECTRUM_COLOR = :red
const BOUNDARY_COLOR = :black
stable_cgrad() = cgrad([STABLE_COLOR, STABLE_COLOR])

# Plot the PI stability region around the slow and fast clusters.
function plot_PI_slow_fast(; cfg, S,
    outdir = joinpath(@__DIR__, "plots", "stability_regions"),
    filename = nothing,
    rads = (slow = nothing, fast = nothing),
    rad_factor = 1.3,
    n = 600,
    stride = 12,
    which = :both,
    spectrum_mode = :unscaled)

    if !(spectrum_mode in (:scaled, :unscaled))
        error("spectrum_mode must be :scaled or :unscaled")
    end

    λ_inner = cfg.dt_inner / cfg.dt

    _, A_prk, b_prk = create_PRK_tableau(
        cfg.c,
        cfg.A,
        cfg.b,
        cfg.K,
        λ_inner
    )

    μ_slow, μ_fast = split_slow_fast(S)

    if spectrum_mode == :scaled
        spec_slow = cfg.dt .* μ_slow
        spec_fast = cfg.dt .* μ_fast

        eval_factor = 1.0

        xlabel_text = "Re(z)"
        ylabel_text = "Im(z)"
        mode_label = "scaled"
    else

        spec_slow = μ_slow
        spec_fast = μ_fast

        eval_factor = cfg.dt

        xlabel_text = "Re(λ)"
        ylabel_text = "Im(λ)"
        mode_label = "unscaled"
    end

    clusters = [spec_slow, spec_fast]
    titles   = ["Slow cluster", "Fast cluster"]
    centers  = mean.(clusters)

    auto_rads = [
        rad_factor * max(maximum(abs.(clusters[j] .- centers[j])), 1.0)
        for j in 1:2
    ]

    plot_rads = [
        rads.slow === nothing ? auto_rads[1] : rads.slow,
        rads.fast === nothing ? auto_rads[2] : rads.fast
    ]

    if which == :slow
        idxs = [1]
    elseif which == :fast
        idxs = [2]
    elseif which == :both
        idxs = [2, 1]
    else
        error("which must be :both, :slow, or :fast")
    end

    plots = Plots.Plot[]

    for j in idxs
        spec = clusters[j]
        center = centers[j]
        rad = plot_rads[j]

        spec_plot = length(spec) > stride ? spec[1:stride:end] : spec

        xs = range(real(center) - rad, real(center) + rad, length = n)
        ys = range(imag(center) - rad, imag(center) + rad, length = n)

        Z_plot = [x + im * y for y in ys, x in xs]

        Z_eval = eval_factor .* Z_plot

        R = abs.(R_grid(A_prk, b_prk, Z_eval))
        stable_fill = ifelse.(R .<= 1.0, 1.0, NaN)

        plt = Plots.heatmap(
            xs, ys, stable_fill;
            color = stable_cgrad(),
            colorbar = false,
            aspect_ratio = :equal,
            xlabel = xlabel_text,
            ylabel = ylabel_text,
            title = titles[j],
            guidefontsize = 11,
            tickfontsize = 9,
            legend = false,
            grid = true,
            framestyle = :box
        )

        Plots.contour!(
            plt, xs, ys, R;
            levels = [1.0],
            color = BOUNDARY_COLOR,
            linewidth = 2,
            label = false
        )

        Plots.scatter!(
            plt, real.(spec_plot), imag.(spec_plot);
            color = SPECTRUM_COLOR,
            markersize = 3,
            markerstrokewidth = 0,
            label = false
        )

        Plots.xlims!(plt, real(center) - rad, real(center) + rad)
        Plots.ylims!(plt, imag(center) - rad, imag(center) + rad)

        push!(plots, plt)
    end

    plt_all = Plots.plot(
        plots...;
        layout = length(plots) == 1 ? (1, 1) : (1, 2),
        size = length(plots) == 1 ? (400, 360) : (860, 360),
        margin = 2Plots.mm,
        bottom_margin = 8Plots.mm,
        left_margin = 5Plots.mm,
        right_margin = 5Plots.mm
    )

    mkpath(outdir)

    if filename === nothing
        dt_label = replace(string(round(cfg.dt, sigdigits = 5)), "." => "p")
        which_label = which == :both ? "slow_fast" : string(which)

        plotfile = joinpath(
            outdir,
            "PI_$(cfg.name)_dt$(dt_label)_$(mode_label)_$(which_label)_clusters.png"
        )
    else
        plotfile = joinpath(outdir, filename)
    end

    savefig(plt_all, plotfile)
    println("Saved plot: ", abspath(plotfile))

    return plt_all
end


# Plot the full PI stability region and the whole spectrum.
function plot_PI_full_spectrum(; cfg, S,
    outdir = joinpath(@__DIR__, "plots", "stability_regions"),
    filename = nothing,
    rad = nothing,
    rad_factor = 1.3,
    n = 600,
    stride = 12,
    spectrum_mode = :unscaled)

    if !(spectrum_mode in (:scaled, :unscaled))
        error("spectrum_mode must be :scaled or :unscaled")
    end

    λ_inner = cfg.dt_inner / cfg.dt

    _, A_prk, b_prk = create_PRK_tableau(
        cfg.c,
        cfg.A,
        cfg.b,
        cfg.K,
        λ_inner
    )

    if spectrum_mode == :scaled
        spec = cfg.dt .* ComplexF64.(S)

        eval_factor = 1.0
        xlabel_text = "Re(z)"
        ylabel_text = "Im(z)"
        mode_label = "scaled"
    else
        spec = ComplexF64.(S)

        eval_factor = cfg.dt
        xlabel_text = "Re(λ)"
        ylabel_text = "Im(λ)"
        mode_label = "unscaled"
    end

    center = mean(spec)

    auto_rad = rad_factor * max(
        maximum(abs.(spec .- center)),
        1.0
    )

    plot_rad = rad === nothing ? auto_rad : rad

    spec_plot = length(spec) > stride ? spec[1:stride:end] : spec

    xs = range(real(center) - plot_rad, real(center) + plot_rad, length = n)
    ys = range(imag(center) - plot_rad, imag(center) + plot_rad, length = n)

    Z_plot = [x + im * y for y in ys, x in xs]
    Z_eval = eval_factor .* Z_plot

    R = abs.(R_grid(A_prk, b_prk, Z_eval))
    stable_fill = ifelse.(R .<= 1.0, 1.0, NaN)

    plt = Plots.heatmap(
        xs, ys, stable_fill;
        color = stable_cgrad(),
        colorbar = false,
        aspect_ratio = :equal,
        xlabel = xlabel_text,
        ylabel = ylabel_text,
        title = "$(cfg.name): full spectrum",
        guidefontsize = 11,
        tickfontsize = 9,
        legend = false,
        grid = true,
        framestyle = :box
    )

    Plots.contour!(
        plt, xs, ys, R;
        levels = [1.0],
        color = BOUNDARY_COLOR,
        linewidth = 2,
        label = false
    )

    Plots.scatter!(
        plt, real.(spec_plot), imag.(spec_plot);
        color = SPECTRUM_COLOR,
        markersize = 3,
        markerstrokewidth = 0,
        label = false
    )
    Plots.xlims!(plt, real(center) - plot_rad, real(center) + plot_rad)
    Plots.ylims!(plt, imag(center) - plot_rad, imag(center) + plot_rad)
    
    mkpath(outdir)

    if filename === nothing
        dt_label = replace(string(round(cfg.dt, sigdigits = 5)), "." => "p")

        plotfile = joinpath(
            outdir,
            "PI_$(cfg.name)_dt$(dt_label)_$(mode_label)_full_spectrum.png"
        )
    else
        plotfile = joinpath(outdir, filename)
    end

    savefig(plt, plotfile)
    println("Saved plot: ", abspath(plotfile))

    return plt
end


# Create and save a shared legend for the PI stability-region plots.
function save_PI_legend(;
    outdir = joinpath(@__DIR__, "plots", "stability_regions"),
    filename = "PI_legend.pdf",
    fontsize = 10,
    ncols = 3,
    size = (760, 55),
    spectrum_mode = :unscaled)

    if !(spectrum_mode in (:scaled, :unscaled))
        error("spectrum_mode must be :scaled or :unscaled")
    end

    if spectrum_mode == :scaled
        spectrum_label = "Scaled spectrum"
        boundary_label = "|R(z)| = 1"
    else
        spectrum_label = "Spectrum"
        boundary_label = "|R(Δt λ)| = 1"
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
        legend_columns = ncols,
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
        label = "PI stable region"
    )

    Plots.plot!(
        leg, [NaN], [NaN];
        color = BOUNDARY_COLOR,
        linewidth = 2,
        label = boundary_label
    )

    Plots.scatter!(
        leg, [NaN], [NaN];
        color = SPECTRUM_COLOR,
        markersize = 4,
        markerstrokewidth = 0,
        label = spectrum_label
    )

    mkpath(outdir)
    legendfile = joinpath(outdir, filename)

    savefig(leg, legendfile)
    println("Saved legend: ", abspath(legendfile))

    return leg
end


