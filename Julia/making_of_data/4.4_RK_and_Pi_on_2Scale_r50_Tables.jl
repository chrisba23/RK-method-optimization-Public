# Reproduces the data for the two-scale problem with r=50 in Chapter 4.4
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

# Two-scale cluster spectrum with r_fast = 50
S = make_eig_clusters(
    α = α,
    ε = ε,
    kappa_slow = 0.05,
    kappa_fast = 0.05
)


# ============================================================
# Table 4.23: Maximal stable step sizes of PFE and PRK4
# ============================================================

println()
println("Table 4.23")
for K in 2:4
    println()
    println("K = $K.")

    outer_tableaus = [
        (
            name = "PFE",
            c = cout_PFE,
            A = Aout_PFE,
            b = bout_PFE,
            K = K,
            dt_inner = dt_inner,
        ),
        (
            name = "PRK4",
            c = cout_RK4,
            A = Aout_RK4,
            b = bout_RK4,
            K = K,
            dt_inner = dt_inner,
        ),
    ]

    compare_PI_methods_on_cluster(
        outer_tableaus = outer_tableaus,
        S = S,
        dt_max_search = 5.0,
    )
end


# ============================================================
# Table 4.24: Optimized 3S* Methods with chebyshev Basis at Polyopt stepsize
# ============================================================

matdir_3S_chebyshev = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebychev", "2-Scale,r=50")
algs_3S_chebyshev = load_algs_auto(matdir_3S_chebyshev; method_type = "3Sstar");
dt_map_3S_chebyshev = load_dt_map_from_polyopt_data(matdir_3S_chebyshev)

println()
println("Table 4.24: Optimized 3S* methods with Chebyshev basis at Polyopt stepsize.")
data_stab = compare_stability_on_spectrum(algs_3S_chebyshev;
    S = S,
    dt_map = dt_map_3S_chebyshev,
    stab_tol = 1e-8
);

# ============================================================
# Table 4.25: Local search around the PolyOpt step size of selected methods.
# ============================================================
println()
println("Table 4.25: Local grid search for p = 1, s = 7.")
diag = diagnose_local_stability(algs_3S_chebyshev;
    S = S,
    dt_map_poly = dt_map_3S_chebyshev,
    method_file = "3Sstar-1-7_2026-04-03T14-02-20.txt",
    factor_left = 0.999,
    factor_right = 1.005,
    n_factors = 30001,
    stab_tol = 0
)

println()
println("Table 4.25: Local grid search for p = 3, s = 9.")
diag = diagnose_local_stability(algs_3S_chebyshev;
    S = S,
    dt_map_poly = dt_map_3S_chebyshev,
    method_file = "3Sstar-3-9_2026-04-22T22-39-05.txt",
    factor_left = 0.999,
    factor_right = 1.005,
    n_factors = 30001,
    stab_tol = 0
)


# ============================================================
#  Optimized 3S* Methods with chebyshev Basis and monotone c at Polyopt stepsize
# ============================================================
println()
println("Optimized 3S* methods with Chebyshev basis and monotone c at Polyopt stepsize.")
matdir_3S_chebyshev_mono = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebychev", "2-Scale,r=50,monotone")
algs_3S_chebyshev_mono = load_algs_auto(matdir_3S_chebyshev_mono; method_type = "3Sstar");
dt_map_3S_chebyshev_mono = load_dt_map_from_polyopt_data(matdir_3S_chebyshev_mono)

data_stab = compare_stability_on_spectrum(algs_3S_chebyshev_mono;
    S = S,
    dt_map = dt_map_3S_chebyshev_mono,
    stab_tol = 1e-8
);

# ============================================================
# Values Reported below Equation (4.68):
# Comparison of Stability-Function Evaluations
# ============================================================

# Compute the maximal stability values for the three evaluations.
function compare_stability_evaluations(alg, zvals; method_type)

    A = alg.coefficients.A
    b = alg.coefficients.b

    R_ls = maximum(abs, R_lowstorage(alg, method_type, zvals))
    R_b  = maximum(abs, R_grid(A, b, zvals))

    coeffs = stability_poly_coeffs(A, b)
    R_p = maximum(z -> abs(evalpoly(z, coeffs)), zvals)

    return (
        low_storage = R_ls,
        butcher = R_b,
        power_series = R_p,
    )
end

matdir_3S_chebyshev = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebychev", "2-Scale,r=50")
method_file = "3Sstar-2-8_2026-05-12T13-18-55.txt"
alg = load_method(matdir_3S_chebyshev, method_file, "3Sstar")
experiment = load_experiment(matdir_3S_chebyshev, method_file)

# Scale the cluster with the PolyOpt step size.
h = Float64(experiment["h_star"])
S_fast = last(split_slow_fast(S))

z_fast = h .* S_fast

# Compare the three stability-function evaluations on the scaled cluster.
stability_evaluation_data = compare_stability_evaluations(
    alg,
    z_fast;
    method_type = "3Sstar"
)
println()
println("Values reported below Equation (4.68) for the scaled fast cluster:")
println("Low-storage evaluation: ", stability_evaluation_data.low_storage)
println("Butcher evaluation:     ", stability_evaluation_data.butcher)
println("Power-series evaluation:", stability_evaluation_data.power_series)

# ============================================================
# Table 4.26: Comparison Float64 and BigFloat evaluations of the stability function
# ============================================================
matdir_3S_chebyshev = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebychev", "2-Scale,r=50")
method_file = "3Sstar-3-10_2026-04-22T22-48-46.txt"

alg = load_method(matdir_3S_chebyshev, method_file, "3Sstar")
experiment = load_experiment(matdir_3S_chebyshev, method_file)

dt = Float64(experiment["h_star"])
S_fast = last(split_slow_fast(S))

z_fast = dt .* S_fast

println()
println("Table 4.26: Comparison of Float64 and BigFloat evaluations of the stability function on the scaled cluster.")
float_comparison = compare_3Sstar_float64_bigfloat(alg, z_fast)



# ============================================================
# Table 4.27: Optimized 2N Methods with chebyshev Basis at Polyopt stepsize
# ============================================================

matdir_2N_chebyshev = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_2N", "Test_2N_chebychev", "2-Scale,r=50")
algs_2N_chebyshev = load_algs_auto(matdir_2N_chebyshev; method_type = "2N")
dt_map_2N_chebyshev = load_dt_map_from_polyopt_data(matdir_2N_chebyshev)

println()
println("Table 4.27: Optimized 2N methods with Chebyshev basis at Polyopt stepsize.")
data_stab = compare_stability_on_spectrum(algs_2N_chebyshev;
    S = S,
    dt_map = dt_map_2N_chebyshev,
    stab_tol = 1e-8
);


# ============================================================
# Table 4.28: Local search around the PolyOpt step size of selected methods.
# ============================================================

println()
println("Table 4.28: Local grid search for p = 3, s = 10.")
diag = diagnose_local_stability(algs_2N_chebyshev;
    S = S,
    dt_map_poly = dt_map_2N_chebyshev,
    method_file = "2N-3-10_2026-04-20T21-53-45.txt",
    factor_left = 0.99999,
    factor_right = 1.000,
    n_factors = 30001,
    stab_tol = 1e-8
)




# ============================================================
# Table 4.29: Optimized 2N Methods with chebyshev Basis and monotone c at Polyopt stepsize
# ============================================================

matdir_2N_chebyshev_mono = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_2N", "Test_2N_chebychev_monotone", "2-Scale,r=50")
algs_2N_chebyshev_mono = load_algs_auto(matdir_2N_chebyshev_mono; method_type = "2N")
dt_map_2N_chebyshev_mono = load_dt_map_from_polyopt_data(matdir_2N_chebyshev_mono)

println()
println("Table 4.29: Optimized 2N methods with Chebyshev basis and monotone c at Polyopt stepsize.")
data_stab = compare_stability_on_spectrum(algs_2N_chebyshev_mono;
    S = S,
    dt_map = dt_map_2N_chebyshev_mono,
    stab_tol = 1e-8
);



# ============================================================
# Table 4.30: Optimized 3S* methods with Interpolation Basis located near the fast cluster
# ============================================================

matdir_interp_fast = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=50","3Sstar_fast_focused")
algs_interp_fast = load_algs_auto(matdir_interp_fast; method_type = "3Sstar")
dt_map_interp_fast = load_dt_map_from_polyopt_data(matdir_interp_fast)

println()
println("Table 4.30: Optimized 3S* methods with Interpolation Basis located near the fast cluster.")
data_stab = compare_stability_on_spectrum(algs_interp_fast;
    S = S,
    dt_map = dt_map_interp_fast,
    stab_tol = 1e-8
);


# ============================================================
# Table 4.31: Local search around the PolyOpt step size of selected methods.
# ============================================================

println()
println("Table 4.31: Local grid search for p = 4, s = 14.")
diag = diagnose_local_stability(algs_interp_fast;
    S = S,
    dt_map_poly = dt_map_interp_fast,
    method_file = "3Sstar-4-14_2026-02-17T17-31-55.txt",
    factor_left = 0.99999,
    factor_right = 1.000,
    n_factors = 30001,
    stab_tol = 1e-8
)

# ============================================================
# Table 4.32: Optimized 3S* methods with Interpolation Basis and different
# locations of the interpolation nodes
# ============================================================

# Interpolation conditions spread over the full interval 
matdir_interp_full = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=50","3Sstar_full_interval")
algs_interp_full = load_algs_auto(matdir_interp_full; method_type = "3Sstar")
dt_map_interp_full = load_dt_map_from_polyopt_data(matdir_interp_full)

println()
println("Table 4.32: Interpolation conditions spread over the full interval.")
data_stab = compare_stability_on_spectrum(algs_interp_full;
    S = S,
    dt_map = dt_map_interp_full,
    stab_tol = 1e-8
);


# Interpolation conditions divided between the slow and fast clusters
matdir_interp_split = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=50","3Sstar_cluster_split")
algs_interp_split = load_algs_auto(matdir_interp_split; method_type = "3Sstar")
dt_map_interp_split = load_dt_map_from_polyopt_data(matdir_interp_split)

println()
println("Table 4.32: Interpolation conditions divided between the slow and fast clusters.")
data_stab = compare_stability_on_spectrum(algs_interp_split;
    S = S,
    dt_map = dt_map_interp_split,
    stab_tol = 1e-8
);


# Interpolation conditions placed only on the fast scale
matdir_interp_fast_only = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=50","3Sstar_only_fast")
algs_interp_fast_only = load_algs_auto(matdir_interp_fast_only; method_type = "3Sstar")
dt_map_interp_fast_only = load_dt_map_from_polyopt_data(matdir_interp_fast_only)

println()
println("Table 4.32: Interpolation conditions placed only on the fast scale.")
data_stab = compare_stability_on_spectrum(algs_interp_fast_only;
    S = S,
    dt_map = dt_map_interp_fast_only,
    stab_tol = 1e-8
);




# ============================================================
# Tables 4.33 and 4.34: Repeated RK-Opt runs
# Each data set contains ten RK-Opt runs for the tested combinations
# ============================================================

# Chebyshev basis with 3S* formulation
matdir_stats_cheby = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebychev", "2-Scale,r=50 - Statistic")
algs_stats_cheby = load_algs_auto(matdir_stats_cheby; method_type = "3Sstar")
dt_map_stats_cheby = load_dt_map_from_polyopt_data(matdir_stats_cheby)

println()
println("Tables 4.33 and 4.34: Repeated RK-Opt runs using the Chebyshev basis with 3S*.")
data_stab = compare_stability_on_spectrum(algs_stats_cheby;
    S = S,
    dt_map = dt_map_stats_cheby,
    stab_tol = 1e-8
);

# Chebyshev interpolation with 3S* formulation fast-focused
matdir_stats_interp = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=50","2-Scale,r=50 - Statistic")
algs_stats_interp = load_algs_auto(matdir_stats_interp; method_type = "3Sstar")
dt_map_stats_interp = load_dt_map_from_polyopt_data(matdir_stats_interp)

println()
println("Tables 4.33 and 4.34: Repeated RK-Opt runs using Chebyshev interpolation with 3S* fast-focused.")
data_stab = compare_stability_on_spectrum(algs_stats_interp;
    S = S,
    dt_map = dt_map_stats_interp,
    stab_tol = 1e-8
);



# Chebyshev basis with 2N formulation
matdir_stats_2N = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_2N", "Test_2N_chebychev", "2-Scale,r=50 - Statistic")
algs_stats_2N = load_algs_auto(matdir_stats_2N; method_type = "2N")
dt_map_stats_2N = load_dt_map_from_polyopt_data(matdir_stats_2N)

println()
println("Tables 4.33 and 4.34: Repeated RK-Opt runs using the Chebyshev basis with 2N.")
data_stab = compare_stability_on_spectrum(algs_stats_2N;
    S = S,
    dt_map = dt_map_stats_2N,
    stab_tol = 1e-8
);


# Chebyshev basis with 2N formulation and monotone c
matdir_stats_2N_monotone = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_2N", "Test_2N_chebychev_monotone", "2-Scale,r=50 - Statistic")
algs_stats_2N_monotone = load_algs_auto(matdir_stats_2N_monotone; method_type = "2N")
dt_map_stats_2N_monotone = load_dt_map_from_polyopt_data(matdir_stats_2N_monotone)

println()
println("Tables 4.33 and 4.34: Repeated RK-Opt runs using the Chebyshev basis with 2N and monotone c.")

data_stab = compare_stability_on_spectrum(algs_stats_2N_monotone;
    S = S,
    dt_map = dt_map_stats_2N_monotone,
    stab_tol = 1e-8
);








# ============================================================
# Table 4.35: Comparison of Butcher and low-storage evaluations
# ============================================================

println()
println("Table 4.35: Comparison of Butcher and low-storage stability evaluations.")
selected_methods = [
    "3Sstar-3-10_2026-04-22T15-36-12.txt",
    "3Sstar-3-11_2026-04-22T19-02-03.txt",
    "3Sstar-4-13_2026-07-23T10-30-49.txt",
    "3Sstar-2-7_2026-04-26T13-00-16.txt",
    "3Sstar-1-8_2026-06-30T11-58-44.txt",
]

print_stability_evaluation_table(
    matdir_interp_fast,
    selected_methods,
    S,
);



# ============================================================
# Table 4.36: Optimized 3S* methods using Chebyshev interpolation
# with interpolation points closest to the fast cluster and the
# low-storage evaluation in the RK-Opt matching conditions
# ============================================================

println()
println("Optimized 3S* methods using Chebyshev interpolation with interpolation points closest to the fast cluster and the low-storage evaluation in the RK-Opt matching conditions.")
matdir_interp_ls  = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=50","3Sstar_fast_focused_low_storage_evaluation")
algs_interp_ls = load_algs_auto(matdir_interp_ls; method_type = "3Sstar")
dt_map_interp_ls = load_dt_map_from_polyopt_data(matdir_interp_ls)

data_stab = compare_stability_on_spectrum(algs_interp_ls;
    S = S,
    dt_map = dt_map_interp_ls,
    stab_tol = 1e-8
);

println()
println("Table 4.36: Comparison of Butcher and low-storage evaluations for the selected methods.")
selected_methods = [
    "3Sstar-3-10_2026-07-04T21-10-56.txt",
    "3Sstar-3-11_2026-07-05T09-57-58.txt",
    "3Sstar-4-11_2026-07-05T12-21-34.txt",
    "3Sstar-4-12_2026-07-19T18-33-19.txt",
]

print_stability_evaluation_table(
    matdir_interp_ls,
    selected_methods,
    S,
);




# ============================================================
# Selected Methods for Tables 4.38-4.39
# ============================================================
algs_all = vcat(algs_interp_fast, algs_interp_ls)
dt_map_all = merge(dt_map_interp_fast, dt_map_interp_ls)

selected_files = [
    "3Sstar-1-8_2026-06-30T11-58-44.txt",
    "3Sstar-2-7_2026-04-26T13-00-16.txt",
    "3Sstar-3-11_2026-07-05T09-57-58.txt",
    "3Sstar-4-13_2026-07-23T10-30-49.txt",
]

algs_selected, dt_map_selected = select_methods(
    algs_all,
    dt_map_all,
    selected_files,
);

println()
println("Table 4.37: Selected 3S* methods for comparison with PI methods.")
data_stab = compare_stability_on_spectrum(
    algs_selected;
    S = S,
    dt_map = dt_map_selected,
    stab_tol = 1e-8
);


# Compute the maximal stable outer step sizes by bisection
K = 3
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
# Table 4.38: Error Comparison for T = 22
# ============================================================

T = 22
println()
println("Table 4.38, RK-Opt methods: Errors for T = 22.")
data_ode_1 = compare_ode(algs_selected;
    problem = :cluster,
    S = S,
    dt_map = dt_map_selected,
    α = α,
    ε = ε,
    tspan = (0.0, T),
    stab_tol = 1e-8
)

println()
println("Table 4.38, PI methods: Errors for T = 22.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 3, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 3, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    problem = :cluster,
    α = α,
    ε = ε,
    S = S,
    tspan = (0.0, T)
)



# ============================================================
# Table 4.39: Error Comparison for T = 100
# ============================================================

T = 100
println()
println("Table 4.39, RK-Opt methods: Errors for T = 100.")
data_ode_1 = compare_ode(algs_selected;
    problem = :cluster,
    S = S,
    dt_map = dt_map_selected,
    α = α,
    ε = ε,
    tspan = (0.0, T),
    stab_tol = 1e-8
)

println()
println("Table 4.39, PI methods: Errors for T = 100.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 3, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 3, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    problem = :cluster,
    α = α,
    ε = ε,
    S = S,
    tspan = (0.0, T)
)

# ============================================================
# Table 4.40: PI Methods with Reduced Step Sizes
# ============================================================
dt_PFE  = 1.15
dt_PRK4 = 1.07
T = 100
println()
println("Table 4.40: PI Methods with Reduced Step Sizes.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 3, dt_inner = dt_inner, dt = dt_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 3, dt_inner = dt_inner, dt = dt_PRK4, p = 4),
    ],
    problem = :cluster,
    α = α,
    ε = ε,
    S = S,
    tspan = (0.0, T)
)

