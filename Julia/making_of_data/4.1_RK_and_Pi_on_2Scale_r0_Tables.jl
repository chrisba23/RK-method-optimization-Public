# Reproduce tables and numerical values from Section 4.1
include(joinpath(@__DIR__, "..", "RK_Data_Functions.jl"))
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

matdir = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=0", "3Sstar_fast_focused")
method_type = "3Sstar"
algs = load_algs_auto(matdir; method_type = method_type)
dt_map_poly = load_dt_map_from_polyopt_data(matdir)


# ============================================================
# Table 4.1: Maximal stable step sizes of PFE and PRK4
# ============================================================

println()
println("Table 4.1: Maximal stable outer step sizes of PFE and PRK4.")
stab_data = compare_PI_methods_on_cluster(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = K, dt_inner = dt_inner),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = K, dt_inner = dt_inner),
    ],
    S = S,
    dt_max_search = 5.0
)


# ============================================================
# Local PolyOpt Stability Diagnosis
# ============================================================

# Evaluate max|R(h λ)| for the PolyOpt polynomial.
function polyopt_stability_at_step(experiment; S, h)

    return maximum(abs.(eval_polyopt_reference(experiment, h .* S)))
end


# Find the connected stable interval containing h_star.
function polyopt_stable_interval(experiment; S,
    rel_window = 2e-10,
    n = 20_001,
    stab_tol = 0.0)


    hstar = Float64(experiment["h_star"])

    factors = range(1 - rel_window, 1 + rel_window; length = n)

    stable = [
        polyopt_stability_at_step(
            experiment;
            S = S,
            h = hstar * factor
        ) <= 1 + stab_tol
        for factor in factors
    ]

    i0 = argmin(abs.(factors .- 1.0))
    stable[i0] || error("hstar is not stable.")

    i_left = i0
    while i_left > 1 && stable[i_left - 1]
        i_left -= 1
    end

    i_right = i0
    while i_right < n && stable[i_right + 1]
        i_right += 1
    end

    i_left == 1 && error("Increase rel_window.")
    i_right == n && error("Increase rel_window.")

    hminus = hstar * factors[i_left]
    hplus  = hstar * factors[i_right]

    return (
        hminus = hminus,
        hplus = hplus,
        delta_minus = hstar - hminus,
        delta_plus = hplus - hstar
    )
end


# Evalutes the local Lipschitz constant of F on the stable interval.
function local_polyopt_lipschitz(experiment; S,
    rel_window = 2e-10,
    n_interval = 20_000,
    n = 10_000,
    stab_tol = 0.0)

    interval = polyopt_stable_interval(
        experiment;
        S = S,
        rel_window = rel_window,
        n = n_interval,
        stab_tol = stab_tol
    )

    hs = range(interval.hminus, interval.hplus; length = n)

    values = [
        polyopt_stability_at_step(
            experiment;
            S = S,
            h = h
        )
        for h in hs
    ]

    return maximum(abs.(diff(values) ./ diff(hs)))
end


# ============================================================
# Equations (4.19)-(4.24): Local Stability near h_star
# ============================================================
method_file = "3Sstar-1-3_2026-06-09T19-43-26.txt"
experiment = load_experiment(matdir, method_file)
h_star = Float64(experiment["h_star"])

println()
println("Equations (4.19)-(4.24)")
# Stability values at the PolyOpt step size.
R_at_hstar = abs.(eval_polyopt_reference(experiment, h_star .* S))
println("|R(-h)| = $(R_at_hstar[1]),   |R(-1000h)| = $(R_at_hstar[2])")

# Test a slightly larger step size.
k = h_star * 1e-10
h_test = h_star + k

println("k = ", k)
# Evaluate the PolyOpt polynomial at the increased step size h + k.
maxR_test = polyopt_stability_at_step(
    experiment;
    S = S,
    h = h_test
)

println(
    "max_{λ ∈ Λ₀} |R((h + k)λ)| = ",
    maxR_test
)

# Sampled stable interval around hstar.
interval = polyopt_stable_interval(
    experiment;
    S = S
)

println(
    "I_loc = [h - $(interval.delta_minus), h + $(interval.delta_plus)]"
)

# Estimate the local Lipschitz constant of F on the stable interval.
L_estimate = local_polyopt_lipschitz(
    experiment;
    S = S
)

println("Local Lipschitz estimate: ", L_estimate)

# ============================================================
# Table 4.7: Optimized 3S* Methods
# ============================================================
println()
println("Table 4.7: Optimized 3S* methods at the PolyOpt step sizes.")

data_stab = compare_stability_on_spectrum(algs;
    S = S,
    dt_map = dt_map_poly,
    stab_tol = 1e-8
);

# Stepsize of p = 1, s = 3 method
method_file = "3Sstar-1-3_2026-06-09T19-43-26.txt"
dt = dt_map_poly[method_file]
println()
println("Stepsize of p = 1, s = 3 method")
println("dt = ", dt)


# Stepsize of p = 1, s = 4 method
method_file = "3Sstar-1-4_2026-07-14T20-55-26.txt"
dt = dt_map_poly[method_file]
println()
println("Stepsize of p = 1, s = 4 method")
println("dt = ", dt)

# Stepsize of p = 2, s = 5 method
method_file = "3Sstar-2-5_2026-06-09T11-07-40.txt"
dt = dt_map_poly[method_file]
println()
println("Stepsize of p = 2, s = 5 method")
println("dt = ", dt)


# ============================================================
# Selected Methods for Tables 4.8--4.10 used in the error comparisons.
# ============================================================
selected_files = [
    "3Sstar-1-4_2026-07-14T20-55-26.txt", 
    "3Sstar-2-4_2026-07-21T19-11-36.txt", 
    "3Sstar-3-4_2026-06-09T18-30-31.txt", 
    "3Sstar-4-5_2026-06-09T18-50-42.txt", 
]

algs_selected, dt_map_selected = select_methods(
    algs,
    dt_map_poly,
    selected_files
);
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


# ============================================================
# Table 4.8: Error Comparison for T = 2.8
# ============================================================
 T = 2.8

println()
println("Table 4.8, RK-Opt methods: Errors for T = 2.8.")
data_ode_1 = compare_ode(algs_selected;
    problem = :twoscale,
    S = S,
    dt_map = dt_map_selected,
    α = α,
    ε = ε,
    tspan = (0.0, T),
    stab_tol = 1e-8
);

println()
println("Table 4.8, PI methods: Errors for T = 2.8.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 1, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 1, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    problem = :twoscale,
    α = α,
    ε = ε,
    tspan = (0.0, T)
);
# ============================================================
# Table 4.9: Error Comparison for T = 100
# ============================================================
println()
println("Table 4.9, RK-Opt methods: Errors for T = 100.")
T = 100
data_ode_1 = compare_ode(algs_selected;
    problem = :twoscale,
    S = S,
    dt_map = dt_map_selected,
    α = α,
    ε = ε,
    tspan = (0.0, T),
    stab_tol = 1e-8
);

println()
println("Table 4.9, PI methods: Errors for T = 100.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 1, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 1, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    problem = :twoscale,
    α = α,
    ε = ε,
    tspan = (0.0, T)
);

# ============================================================
# Table 4.10: Error Comparison for T = 32128
# ============================================================

T = 32128

println()
println("Table 4.10, RK-Opt methods: Errors for T = 32128.")
data_ode_2 = compare_ode(algs_selected;
    problem = :twoscale,
    S = S,
    dt_map = dt_map_selected,
    α = α,
    ε = ε,
    tspan = (0.0, T),
    stab_tol = 1e-8
);

println()
println("Table 4.10, PI methods: Errors for T = 32128.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 1, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 1, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    problem = :twoscale,
    α = α,
    ε = ε,
    tspan = (0.0, T)
);


# ============================================================
# Table 4.11: PI Methods with Reduced Step Sizes
# ============================================================
dt_PFE  = 1.2
dt_PRK4 = 1.12
T = 100
println()
println("Table 4.11: Errors of PFE and PRK4 with reduced outer step sizes.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 1, dt_inner = dt_inner, dt = dt_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 1, dt_inner = dt_inner, dt = dt_PRK4, p = 4),
    ],
    problem = :twoscale,
    α = α,
    ε = ε,
    tspan = (0.0, T)
);

