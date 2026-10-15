# ============================================================
# Projective integration stability and ODE experiments
#
# Defines the outer Runge-Kutta methods used in the
# projective schemes, computes stability limits on a given
# spectrum, and compares the methods on the ODE test problems.
# ============================================================

include(joinpath(@__DIR__, "setup.jl"))

# Check PI stability for a given outer step size.
function stable_for_dt_cluster_PI(dt; K, dt_inner, S, c_out, A_out, b_out,
                                  stab_tol = 1e-8)

    dt <= (K + 1) * dt_inner && return false, Inf

    λ = dt_inner / dt
    _, A_prk, b_prk = create_PRK_tableau(c_out, A_out, b_out, K, λ)

    mx = 0.0

    for μ in S
        mx = max(mx, abs(R_general(A_prk, b_prk, dt * μ)))
    end

    return (mx <= 1.0 + stab_tol), mx
end

# Find the largest stable outer step size on the given spectrum by bisection.
function find_dt_max_cluster_PI(; K, dt_inner, S, c_out, A_out, b_out,
                                dt_max_search = 5.0,
                                stab_tol = 1e-8,
                                bisect_tol = 1e-10,
                                max_iter = 100)

    dt_min = (K + 1) * dt_inner * 1.001 

    lo = dt_min
    hi = dt_max_search

    ok_lo, _ = stable_for_dt_cluster_PI(lo;
        K = K,
        dt_inner = dt_inner,
        S = S,
        c_out = c_out,
        A_out = A_out,
        b_out = b_out,
        stab_tol = stab_tol
    )

    ok_hi, _ = stable_for_dt_cluster_PI(hi;
        K = K,
        dt_inner = dt_inner,
        S = S,
        c_out = c_out,
        A_out = A_out,
        b_out = b_out,
        stab_tol = stab_tol
    )

    ok_lo || error("Lower bound dt_min = $(dt_min) is not stable.")
    ok_hi && error("Upper bound dt_max_search = $(dt_max_search) is still stable. Increase dt_max_search.")

    for _ in 1:max_iter
        mid = 0.5 * (lo + hi)

        ok_mid, _ = stable_for_dt_cluster_PI(mid;
            K = K,
            dt_inner = dt_inner,
            S = S,
            c_out = c_out,
            A_out = A_out,
            b_out = b_out,
            stab_tol = stab_tol
        )

        if ok_mid
            lo = mid
        else
            hi = mid
        end

        if hi - lo < bisect_tol
            break
        end
    end

    _, mx = stable_for_dt_cluster_PI(lo;
        K = K,
        dt_inner = dt_inner,
        S = S,
        c_out = c_out,
        A_out = A_out,
        b_out = b_out,
        stab_tol = stab_tol
    )

    return lo, mx
end


# Compare the stability limits and dt/#RHS of PI methods on the given spectrum..
function compare_PI_methods_on_cluster(; outer_tableaus, S,
                                       dt_max_search = 5.0,
                                       stab_tol = 1e-8,
                                       bisect_tol = 1e-10,
                                       max_iter = 100)

    rows = Any[]

    for cfg in outer_tableaus
        c = cfg.c
        A = cfg.A
        b = cfg.b
        K = cfg.K
        dt_inner = cfg.dt_inner
        name = cfg.name

        dtmax, mx = find_dt_max_cluster_PI(
            K = K,
            dt_inner = dt_inner,
            S = S,
            c_out = c,
            A_out = A,
            b_out = b,
            dt_max_search = dt_max_search,
            stab_tol = stab_tol,
            bisect_tol = bisect_tol,
            max_iter = max_iter)

        rhs_per_step = (K + 1) * length(b)
        dtrhs = dtmax / rhs_per_step

        push!(rows, Any[
            name,
            K,
            dt_inner,
            rhs_per_step,
            dtmax,
            dtrhs,
            mx
        ])
    end

    data = reduce(vcat, permutedims.(rows))

    pretty_table(
        data;
        header = [
            "Algorithm",
            "K",
            "dt_inner",
            "#RHS/step",
            "dt_max",
            "dt/#RHS",
            "max|R|"
        ],
        crop = :none,
        formatters = (
            ft_printf("%s",   [1]),
            ft_printf("%d",   [2]),
            ft_printf("%.3e", [3]),
            ft_printf("%d",   [4]),
            ft_printf("%.6f", [5]),
            ft_printf("%.6f", [6]),
            ft_printf("%.6e", [7])
        ),
    )

    return data
end


# Compare the numerical errors and dt/#RHS of PI methods
# for the two-scale and clustered-spectrum ODE problems.
function compare_PI_methods_ode(; outer_tableaus,
    problem = :twoscale,
    S = nothing,
    α = 1.0,
    ε = 1e-3,
    tspan = (0.0, 1e3),
    u10 = 1.0,
    u20 = 0.0,
    u0 = nothing,
    stab_tol = 1e-8,
    show_unstable = false)

    rows = Any[]
    T = tspan[end]

    if problem == :twoscale
        S_check = ComplexF64[-α, -1 / ε]

        slow_idx = [1]
        fast_idx = [2]

        u0_cluster = nothing

    elseif problem == :cluster
        S === nothing && error("For problem = :cluster, provide S.")

        S_check = ComplexF64.(S)
        length(S_check) >= 2 || error("S must contain slow and fast eigenvalues.")
        iseven(length(S_check)) || error("S must contain equally many slow and fast eigenvalues.")

        n_slow = div(length(S_check), 2)
        slow_idx = collect(1:n_slow)
        fast_idx = collect((n_slow + 1):length(S_check))

        u0_cluster = u0 === nothing ?
            ones(ComplexF64, length(S_check)) :
            ComplexF64.(u0)

        length(u0_cluster) == length(S_check) ||
            error("length(u0) must match length(S).")

    else
        error("Use problem = :twoscale or :cluster.")
    end

    for cfg in outer_tableaus
        c = cfg.c
        A = cfg.A
        b = cfg.b
        K = cfg.K
        dt_inner = cfg.dt_inner
        dt = cfg.dt
        p = cfg.p
        name = cfg.name

        λ_inner = dt_inner / dt

        c_prk, A_prk, b_prk = create_PRK_tableau(
            c, A, b, K, λ_inner
        )

        tab = DiffEqBase.ExplicitRKTableau(A_prk, c_prk, b_prk, p)
        alg = OrdinaryDiffEq.ExplicitRK(tab)

        ok, mx = stable_for_dt_cluster_PI(dt;
            K = K,
            dt_inner = dt_inner,
            S = S_check,
            c_out = c,
            A_out = A,
            b_out = b,
            stab_tol = stab_tol
        )

        if !show_unstable && !ok
            continue
        end


        if problem == :twoscale
            u1T_ex, u2T_ex = exact_two_scale(T;
                α = α,
                ε = ε,
                u10 = u10,
                u20 = u20
            )

            u_exact = [u1T_ex, u2T_ex]

            sol = solve_two_scale(alg;
                dt = dt,
                α = α,
                ε = ε,
                tspan = (0.0, T),
                u10 = u10,
                u20 = u20
            )

        else
            u_exact = exact_cluster(T, S_check; u0 = u0_cluster)

            sol = solve_cluster_ode(alg;
                S = S_check,
                dt = dt,
                tspan = (0.0, T),
                u0 = u0_cluster
            )
        end

        uT = sol.u[end]
        err_vec = uT .- u_exact

        abs_slow = norm(err_vec[slow_idx], Inf)
        abs_fast = norm(err_vec[fast_idx], Inf)

        abs_error = max(abs_slow, abs_fast)
        rel_error = if iszero(norm(u_exact, Inf))
            NaN
        else
            abs_error / norm(u_exact, Inf)
        end

        rhs_per_step = length(b_prk)
        n_rhs = sol.stats.naccept * rhs_per_step

        push!(rows, Any[
            name,
            K,
            dt_inner,
            dt,
            dt / rhs_per_step,
            T,
            n_rhs,
            abs_slow,
            abs_fast,
            abs_error,
            rel_error,
            mx
        ])
    end

    isempty(rows) && error("No methods left in table. Maybe all methods are unstable.")

    sort!(rows, by = r -> r[5], rev = true)
    data = reduce(vcat, permutedims.(rows))

    pretty_table(
        data;
        header = [
            "Algorithm",
            "K",
            "dt_inner",
            "dt",
            "dt/#RHS",
            "T",
            "#RHS",
            "abs slow",
            "abs fast",
            "abs total",
            "rel error",
            "max|R|"
        ],
        crop = :none,
        formatters = (
            ft_printf("%s",       [1]),
            ft_printf("%d",       [2]),
            ft_printf("%.3e",     [3]),
            ft_printf("%.3e",     [4]),
            ft_printf("%.3e",     [5]),
            ft_printf("%.3e",     [6]),
            ft_printf("%d",       [7]),
            ft_printf("%.3e",     [8]),
            ft_printf("%.3e",     [9]),
            ft_printf("%.3e",     [10]),
            ft_printf("%.3e",     [11]),
            ft_printf("%.6e",     [12])
        ),
    )

    return data
end

# Compare the numerical errors and dt/#RHS of PI methods
# for the generalized two-scale linear system.
function compare_PI_methods_generalized_two_scale_ode(;
    outer_tableaus,
    A_sys,
    S,
    u0,
    tspan = (0.0, 1.0),
    stab_tol = 1e-8,
    show_unstable = false)

    rows = Any[]

    n_slow = div(length(u0), 2)
    slow_idx = 1:n_slow
    fast_idx = (n_slow + 1):length(u0)

    T = tspan[end]

    for cfg in outer_tableaus
        c = cfg.c
        A_out = cfg.A
        b = cfg.b
        K = cfg.K
        dt_inner = cfg.dt_inner
        dt = cfg.dt
        p = cfg.p
        name = cfg.name

        λ_inner = dt_inner / dt

        c_prk, A_prk, b_prk = create_PRK_tableau(
            c, A_out, b, K, λ_inner
        )

        tab = DiffEqBase.ExplicitRKTableau(A_prk, c_prk, b_prk, p)
        alg = OrdinaryDiffEq.ExplicitRK(tab)

        ok, mx = stable_for_dt_cluster_PI(
            dt;
            K = K,
            dt_inner = dt_inner,
            S = S,
            c_out = c,
            A_out = A_out,
            b_out = b,
            stab_tol = stab_tol
        )

        if !show_unstable && !ok
            continue
        end

        sol = solve_generalized_two_scale(alg;
            A_sys = A_sys,
            dt = dt,
            tspan = (0.0, T),
            u0 = u0
        )

        u_exact = exp(T * A_sys) * u0
        err_vec = sol.u[end] .- u_exact

        abs_slow = norm(err_vec[slow_idx], Inf)
        abs_fast = norm(err_vec[fast_idx], Inf)

        abs_error = max(abs_slow, abs_fast)
        rel_error = if iszero(norm(u_exact, Inf))
            NaN
        else
            abs_error / norm(u_exact, Inf)
        end

        rhs_per_step = length(b_prk)
        n_rhs = sol.stats.naccept * rhs_per_step

        push!(rows, Any[
            name,
            K,
            dt_inner,
            dt,
            dt / rhs_per_step,
            T,
            n_rhs,
            abs_slow,
            abs_fast,
            abs_error,
            rel_error,
            mx
        ])
    end

    isempty(rows) &&
        error("No PI methods left in table. Maybe all methods are unstable.")

    sort!(rows, by = r -> r[5], rev = true)
    data = reduce(vcat, permutedims.(rows))

    pretty_table(
        data;
        header = [
            "Algorithm",
            "K",
            "dt_inner",
            "dt",
            "dt/#RHS",
            "T",
            "#RHS",
            "abs slow",
            "abs fast",
            "abs total",
            "rel error",
            "max|R|"
        ],
        crop = :none,
        formatters = (
            ft_printf("%s",       [1]),
            ft_printf("%d",       [2]),
            ft_printf("%.3e",     [3]),
            ft_printf("%.3e",     [4]),
            ft_printf("%.3e",     [5]),
            ft_printf("%.3e",     [6]),
            ft_printf("%d",       [7]),
            ft_printf("%.3e",     [8]),
            ft_printf("%.3e",     [9]),
            ft_printf("%.3e",     [10]),
            ft_printf("%.3e",     [11]),
            ft_printf("%.6e",     [12])
        ),
    )

    return data
end

# Create a work-precision diagram for PFE and PRK4
# on the two-scale, clustered-spectrum, or generalized ODE problem.
function work_precision_PI(;
    problem = :twoscale,
    S =  nothing,
    A_sys = nothing,
    u0 = nothing,
    T = 100.0,
    K = 3,
    dt_inner = 1e-3,
    factors = [1.0, 0.8, 0.6, 0.4, 0.2, 0.1],
    annotate = false,
    outdir = joinpath(@__DIR__, "plots"),
    filename = nothing
)

    α = 1.0
    ε = 1e-3

    problem in (:twoscale, :cluster, :generalized) ||
        error("Use problem = :twoscale, :cluster, or :generalized.")

    S === nothing &&
        error("Provide S.")

    if problem == :generalized
        A_sys === nothing &&
            error("For problem = :generalized, provide A_sys.")

        u0 === nothing &&
            error("For problem = :generalized, provide u0.")
    end

    S = ComplexF64.(S)

    methods = [
        (name = "PFE",  p = 1, c = cout_PFE, A = Aout_PFE, b = bout_PFE),
        (name = "PRK4", p = 4, c = cout_RK4, A = Aout_RK4, b = bout_RK4)
    ]

    plt = Plots.plot(
        ;
        xscale = :log10,
        yscale = :log10,
        xlabel = "RHS evaluations",
        ylabel = "Absolute error",
        legend = :topright,
        guidefontsize = 14,
        tickfontsize = 11,
        legendfontsize = 12,
        left_margin = 4Plots.mm,
        right_margin = 8Plots.mm,
        top_margin = 4Plots.mm,
        bottom_margin = 4Plots.mm,
        size = (900, 600)
    )

    for (method, marker) in zip(methods, [:diamond, :star5])

        dtmax, _ = find_dt_max_cluster_PI(
            K = K,
            dt_inner = dt_inner,
            S = S,
            c_out = method.c,
            A_out = method.A,
            b_out = method.b,
            dt_max_search = 5.0
        )

        rhs = Int[]
        errors = Float64[]
        dts = Float64[]

        for factor in factors
            dt = factor * dtmax
            dt > T && continue
            
            λ = dt_inner / dt

            c_prk, A_prk, b_prk = create_PRK_tableau(
                method.c, method.A, method.b, K, λ
            )

            tab = DiffEqBase.ExplicitRKTableau(
                A_prk, c_prk, b_prk, method.p
            )

            alg = OrdinaryDiffEq.ExplicitRK(tab)

            if problem == :twoscale
                sol = solve_two_scale(
                    alg;
                    dt = dt,
                    α = α,
                    ε = ε,
                    tspan = (0.0, T)
                )

                u1, u2 = exact_two_scale(T; α = α, ε = ε)
                u_exact = [u1, u2]

            elseif problem == :cluster
                u0_cluster = ones(ComplexF64, length(S))

                sol = solve_cluster_ode(
                    alg;
                    S = S,
                    dt = dt,
                    tspan = (0.0, T),
                    u0 = u0_cluster
                )

                u_exact = exact_cluster(T, S; u0 = u0_cluster)

            elseif problem == :generalized
                sol = solve_generalized_two_scale(alg;
                    A_sys = A_sys,
                    dt = dt,
                    tspan = (0.0, T),
                    u0 = u0
                )

                u_exact = exp(T * A_sys) * u0
            end

            err = norm(sol.u[end] .- u_exact, Inf)

            push!(rhs, sol.stats.naccept * length(b_prk))
            push!(errors, err)
            push!(dts, dt)
        end

        perm = sortperm(rhs)

        Plots.plot!(
            plt,
            rhs[perm],
            errors[perm];
            label = method.name,
            marker = marker,
            markersize = 8,
            linewidth = 2,
            linestyle = :dash
        )

        if annotate
            for j in perm
                Plots.annotate!(
                    plt,
                    rhs[j],
                    errors[j],
                    Plots.text(
                        "dt=$(round(dts[j]; sigdigits = 3))",
                        9,
                        :left
                    )
                )
            end
        end
    end

    if filename !== nothing
        mkpath(outdir)

        plotfile = joinpath(outdir, filename)
        savefig(plt, plotfile)

        println("Saved plot: ", abspath(plotfile))
    end

    return plt
end