# ============================================================
# Functions for RK method experiments
# This file loads RK methods and PolyOpt step sizes, checks stability on a
# spectrum, compares ODE errors, and performs a local stability diagnosis around the PolyOpt step size.
# ============================================================

include(joinpath(@__DIR__, "setup.jl"))


# Load all RK method files in a folder using the given reader.
function load_algs_from(matdir, reader)
    files = filter(f -> endswith(lowercase(f), ".txt"), readdir(matdir; join=true))
    return [reader(filename=f) for f in files]
end

function load_algs_auto(matdir; method_type)
    if method_type == "3Sstar"
        return load_algs_from(matdir, RK3SFromFile)
    elseif method_type == "2N"
        return load_algs_from(matdir, RK2NFromFile)
    elseif method_type == "ERK"
        return load_algs_from(matdir, ExplicitRKFromFile)
    else
        error("Unknown method_type: $method_type")
    end
end


# Load the PolyOpt step sizes saved for all optimized RK methods.
function load_dt_map_from_polyopt_data(matdir)
    data_dir = joinpath(matdir, "polyopt_data")

    isdir(data_dir) || error("polyopt_data folder not found: $data_dir")

    mat_files = filter(f -> endswith(lowercase(f), "_plot_data.mat"), readdir(data_dir; join=true))
    isempty(mat_files) && error("No .mat files found in: $data_dir")

    dt_map = Dict{String, Float64}()

    for f in mat_files
        vars = matread(f)

        haskey(vars, "experiment") || error("No variable 'experiment' found in $(basename(f))")

        experiment = vars["experiment"]

        method_file = string(experiment["method_file"])
        dt = Float64(experiment["h_star"])

        dt_map[basename(method_file)] = dt
    end

    return dt_map
end

# Select chosen RK methods for comparison.
function select_methods(algs, dt_map, files)

    alg_by_file = Dict(
        basename(alg.filename) => alg
        for alg in algs
    )

    algs_selected = [alg_by_file[file] for file in files]

    dt_map_selected = Dict{String, Float64}()

    for file in files
        dt_map_selected[file] = dt_map[file]
    end

    return algs_selected, dt_map_selected
end

# Check stability for one method on a spectrum and with a step size.
function stable_for_dt_cluster_alg(dt; alg, S, stab_tol=1e-8)
    A = alg.coefficients.A
    b = alg.coefficients.b

    mx = 0.0

    for λeig in S
        mx = max(mx, abs(R_general(A, b, dt * λeig)))
    end

    return (mx <= 1.0 + stab_tol), mx
end

# Read order p and number of stages s from a filename
function method_sort_key(filename)

    parts = split(basename(filename), "-")

    p = parse(Int, parts[2])
    s = parse(Int, split(parts[3], "_")[1])

    return (p, s)
end

# Compare the stability and dt/#RHS of RK methods on the given spectrum.
function compare_stability_on_spectrum(algs; S, dt_map,
    stab_tol = 1e-8)

    rows = Any[]

    for alg in algs
        name = basename(alg.filename)
        haskey(dt_map, name) || error("No dt found for '$name'.")

        dt = dt_map[name]

        ok, mx = stable_for_dt_cluster_alg(dt;
            alg = alg,
            S = S,
            stab_tol = stab_tol
        )

        status = ok ? "stable" : "not stable"

        s = alg.coefficients.s

        push!(rows, Any[
            name,
            dt,
            dt / s,
            mx,
            status
        ])
    end

    sort!(rows, by = r -> method_sort_key(r[1]))
    data = reduce(vcat, permutedims.(rows))

    pretty_table(
        data;
        header = ["Algorithm", "dt", "dt/#RHS", "max|R|", "Status"],
        crop = :none,
        formatters = (
            ft_printf("%s",     [1]),
            ft_printf("%.3e",   [2]),
            ft_printf("%.3e",   [3]),
            ft_printf("%.6e",   [4]),
            ft_printf("%s",     [5])
        ),
    )

    return data
end

# Compare the numerical errors and dt/#RHS of RK methods 
# for the two-scale and clustered-spectrum ODE problems.
function compare_ode(algs; problem, S, dt_map,
    α = 1.0,
    ε = 1e-3,
    tspan = (0.0, 1e3),
    stab_tol = 1e-8,
    u10 = 1.0,
    u20 = 0.0,
    u0 = ones(ComplexF64, length(S)),
    show_unstable = false)

    rows = Any[]

    T_target = tspan[end]


    if problem == :twoscale
        S_stab = ComplexF64[-α, -1 / ε]

        slow_idx = [1]
        fast_idx = [2]

    elseif problem == :cluster
        S_stab = S
        length(u0) == length(S) || error("length(u0) must match length(S).")
        
        if !iseven(length(S))
            error("S must contain equally many slow and fast eigenvalues.")
        end

        n_slow = div(length(S), 2)

        slow_idx = collect(1:n_slow)
        fast_idx = collect((n_slow + 1):length(S))

        isempty(slow_idx) && error("slow_idx is empty.")
        isempty(fast_idx) && error("fast_idx is empty. Check S.")

    else
        error("Unknown problem: $problem. Use :twoscale or :cluster.")
    end

    for alg in algs
        name = basename(alg.filename)
        haskey(dt_map, name) || error("No dt found for '$name'.")

        dt = dt_map[name]

        # Full steps only: run until T_target is reached or exceeded.
        # This avoids a final remainder step.
        nsteps = ceil(Int, T_target / dt)
        T_used = nsteps * dt
        tspan_used = (0.0, T_used)

        ok, mx = stable_for_dt_cluster_alg(dt;
            alg = alg,
            S = S_stab,
            stab_tol = stab_tol
        )

        if !show_unstable && !ok
            continue
        end

        if problem == :twoscale
            u1T_ex, u2T_ex = exact_two_scale(T_used;
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
                tspan = tspan_used,
                u10 = u10,
                u20 = u20
            )

        elseif problem == :cluster
            u_exact = exact_cluster(T_used, S; u0 = u0)

            sol = solve_cluster_ode(alg;
                S = S,
                dt = dt,
                tspan = tspan_used,
                u0 = u0
            )
        end

        s = alg.coefficients.s
        n_rhs = s * sol.stats.naccept

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

        push!(rows, Any[
            name,
            dt,
            dt / s,
            T_used,
            n_rhs,
            abs_slow,
            abs_fast,
            abs_error,
            rel_error,
            mx
        ])
    end

    isempty(rows) && error("No methods left in table. Maybe all methods are unstable.")

    sort!(rows, by = r -> method_sort_key(r[1]))
    data = reduce(vcat, permutedims.(rows))

    pretty_table(
        data;
        header = [
            "Algorithm",
            "dt",
            "dt/#RHS",
            "T used",
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
            ft_printf("%.3e",     [2]),
            ft_printf("%.3e",     [3]),
            ft_printf("%.3e",     [4]),
            ft_printf("%d",       [5]),
            ft_printf("%.3e",     [6]),
            ft_printf("%.3e",     [7]),
            ft_printf("%.3e",     [8]),
            ft_printf("%.3e",     [9]),
            ft_printf("%.6e",     [10])
        ),
    )

    return data
end

# Search around the PolyOpt step size for a stable step size
# of one selected RK method.
function diagnose_local_stability(algs; S, dt_map_poly, method_file,
    factor_left = 0.999990,
    factor_right = 1.000010,
    n_factors = 30001,
    stab_tol = 1e-8)

    matches = filter(a -> basename(a.filename) == method_file, algs)
    length(matches) == 1 || error("Could not uniquely find method file: $method_file")

    alg = only(matches)
    hstar = dt_map_poly[method_file]

    _, mx_h = stable_for_dt_cluster_alg(hstar;
        alg = alg,
        S = S,
        stab_tol = stab_tol
    )

    factors = range(factor_left, factor_right; length = n_factors)

    best = (
        factor = NaN,
        dt = NaN,
        mx = Inf
    )

    stable_best = nothing

    for fac in factors
        dt = hstar * fac

        ok, mx = stable_for_dt_cluster_alg(dt;
            alg = alg,
            S = S,
            stab_tol = stab_tol
        )

        if mx < best.mx
            best = (
                factor = fac,
                dt = dt,
                mx = mx
            )
        end

        if ok
            candidate = (
                factor = fac,
                dt = dt,
                mx = mx
            )

            if stable_best === nothing || candidate.dt > stable_best.dt
                stable_best = candidate
            end
        end
    end

    println("h_PolyOpt        = ", hstar)
    println("max|R(h λ)|      = ", mx_h)

    if stable_best === nothing
        println("No stable local dt found.")
        println("best local factor = ", best.factor)
        println("best local dt     = ", best.dt)
        println("best local max|R| = ", best.mx)
    else
        println("local stable factor = ", stable_best.factor)
        println("local stable dt     = ", stable_best.dt)
        println("max|R(dt λ)|        = ", stable_best.mx)
    end

    return stable_best, best
end


# Compare the numerical errors and dt/#RHS of RK methods
# for the generalized two-scale linear system.
function compare_generalized_two_scale_ode(
    algs;
    A_sys,
    S,
    u0,
    dt_map,
    tspan = (0.0, 0.5),
    stab_tol = 1e-8,
    show_unstable = false)

    rows = Any[]

    n_slow = div(length(u0), 2)
    slow_idx = 1:n_slow
    fast_idx = (n_slow + 1):length(u0)

    T_target = tspan[end]

    for alg in algs
        name = basename(alg.filename)
        haskey(dt_map, name) || error("No dt found for '$name'.")
        dt = dt_map[name]

        nsteps = ceil(Int, T_target / dt)
        T_used = nsteps * dt

        ok, mx = stable_for_dt_cluster_alg(
            dt;
            alg = alg,
            S = S,
            stab_tol = stab_tol,
        )

        if !show_unstable && !ok
            continue
        end

        sol = solve_generalized_two_scale(alg;
            A_sys = A_sys,
            dt = dt,
            tspan = (0.0, T_used),
            u0 = u0
        )
        u_exact = exp(T_used * A_sys) * u0
        err_vec = sol.u[end] .- u_exact

        abs_slow = norm(err_vec[slow_idx], Inf)
        abs_fast = norm(err_vec[fast_idx], Inf)

        abs_error = max(abs_slow, abs_fast)
        rel_error = if iszero(norm(u_exact, Inf))
            NaN
        else
            abs_error / norm(u_exact, Inf)
        end

        s = alg.coefficients.s
        n_rhs = s * sol.stats.naccept

        push!(rows, Any[
            name,
            dt,
            dt / s,
            T_used,
            n_rhs,
            abs_slow,
            abs_fast,
            abs_error,
            rel_error,
            mx
        ])
    end
    isempty(rows) &&
        error("No methods left in table. Maybe all methods are unstable.")
    sort!(rows, by = r -> method_sort_key(r[1]))
    data = reduce(vcat, permutedims.(rows))

    pretty_table(
        data;
        header = [
            "Algorithm",
            "dt",
            "dt/#RHS",
            "T used",
            "#RHS",
            "abs slow",
            "abs fast",
            "abs total",
            "rel error",
            "max|R|"
        ],
        crop = :none,
        formatters = (
            ft_printf("%s",   [1]),
            ft_printf("%.3e", [2]),
            ft_printf("%.3e", [3]),
            ft_printf("%.3e", [4]),
            ft_printf("%d",   [5]),
            ft_printf("%.3e", [6]),
            ft_printf("%.3e", [7]),
            ft_printf("%.3e", [8]),
            ft_printf("%.3e", [9]),
            ft_printf("%.6e", [10])
        ),
    )

    return data
end
# Compare Float64 and BigFloat stability-function evaluations
# in 3S*, power-series, and Butcher form.
function compare_3Sstar_float64_bigfloat(alg, zvals)

    beta, gamma1, gamma2, gamma3, delta =
        extract_3Sstar_params(alg)

    A = alg.coefficients.A
    b = alg.coefficients.b

    z64 = ComplexF64.(zvals)

    Rls64 = R_3Sstar(
        beta,
        gamma1,
        gamma2,
        gamma3,
        delta,
        z64,
    )

    coeffs64 = stability_poly_coeffs(A, b)
    Rp64 = map(z -> evalpoly(z, coeffs64), z64)
    Rb64 = R_grid(A, b, z64)

    err_lowstorage, err_power, err_butcher =
        setprecision(BigFloat, 256) do

        zbig = Complex{BigFloat}.(zvals)
        Abig = BigFloat.(A)
        bbig = BigFloat.(b)

        Rlsbig = R_3Sstar(
            BigFloat.(beta),
            BigFloat.(gamma1),
            BigFloat.(gamma2),
            BigFloat.(gamma3),
            BigFloat.(delta),
            zbig,
        )

        coeffsbig = stability_poly_coeffs(Abig, bbig)
        Rpbig = map(z -> evalpoly(z, coeffsbig), zbig)
        Rbbig = R_grid(Abig, bbig, zbig)

        err_lowstorage_big = maximum(
            abs.(Complex{BigFloat}.(Rls64) .- Rlsbig) ./
            abs.(Rlsbig)
        )

        err_power_big = maximum(
            abs.(Complex{BigFloat}.(Rp64) .- Rpbig) ./
            abs.(Rpbig)
        )

        err_butcher_big = maximum(
            abs.(Complex{BigFloat}.(Rb64) .- Rbbig) ./
            abs.(Rbbig)
        )

        return (
            err_lowstorage_big,
            err_power_big,
            err_butcher_big,
        )
    end

    rows = Any[
        ["Low-storage", Float64(err_lowstorage)],
        ["Power series", Float64(err_power)],
        ["Butcher", Float64(err_butcher)],
    ]

    data = reduce(vcat, permutedims.(rows))

    pretty_table(
        data;
        header = ["Evaluation", "Maximum relative difference"],
        formatters = (
            ft_printf("%s", [1]),
            ft_printf("%.3e", [2])
        ),
    )

    return data
end

# Compare the maximal stability values on the slow and fast clusters.
function clusterwise_stability_RK(algs; S, dt_map)
    n_slow = length(S) ÷ 2
    S_slow = S[1:n_slow]
    S_fast = S[(n_slow + 1):end]

    rows = Any[]

    for alg in algs
        name = basename(alg.filename)
        dt = dt_map[name]

        A = alg.coefficients.A
        b = alg.coefficients.b
        s = alg.coefficients.s

        R_slow_max = maximum(
            abs(R_general(A, b, dt * λ))
            for λ in S_slow
        )

        R_fast_max = maximum(
            abs(R_general(A, b, dt * λ))
            for λ in S_fast
        )

        push!(rows, Any[
            name,
            dt,
            dt / s,
            R_slow_max,
            R_fast_max,
            max(R_slow_max, R_fast_max),
        ])
    end

    sort!(rows, by = row -> method_sort_key(row[1]))
    data = reduce(vcat, permutedims.(rows))

    pretty_table(
        data;
        header = [
            "Algorithm",
            "dt",
            "dt/#RHS",
            "max |R_slow|",
            "max |R_fast|",
            "max |R|",
        ],
        crop = :none,
        formatters = (
            ft_printf("%s", [1]),
            ft_printf("%.3e", [2]),
            ft_printf("%.3e", [3]),
            ft_printf("%.6e", [4]),
            ft_printf("%.6e", [5]),
            ft_printf("%.6e", [6]),
        ),
    )

    return data
end




# Compare the Butcher and low-storage stability evaluations for selected methods 
function print_stability_evaluation_table(
    matdir,
    selected_methods,
    S;
    method_type = "3Sstar",
)
    rows = Any[]

    for method_file in selected_methods
        experiment = load_experiment(matdir, method_file)
        alg = load_method(matdir, method_file, method_type)

        dt = Float64(experiment["h_star"])
        zvals = dt .* S

        A = alg.coefficients.A
        b = alg.coefficients.b

        max_R_butcher = maximum(abs, R_grid(A, b, zvals))
        max_R_low_storage = maximum(
            abs,
            R_lowstorage(alg, method_type, zvals),
        )

        push!(rows, Any[
            method_file,
            max_R_butcher,
            max_R_low_storage,
        ])
    end

    sort!(rows, by = row -> method_sort_key(row[1]))
    data_table = reduce(vcat, permutedims.(rows))

    pretty_table(
        data_table;
        header = [
            "Algorithm",
            "max |R_B|",
            "max |R_LS|",
        ],
        crop = :none,
        formatters = (
            ft_printf("%s", [1]),
            ft_printf("%.6e", [2]),
            ft_printf("%.6e", [3]),
        ),
    )

    return data_table
end