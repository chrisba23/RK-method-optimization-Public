# Reproduces the data for the generalized two-scale problem in Chapter 4.6
include(joinpath(@__DIR__, "..", "RK_Data_Functions.jl"))
include(joinpath(@__DIR__, "..", "Pi_Methods_Data.jl"))

# ============================================================
# Load Problem Matrix
# ============================================================

matrixdir = joinpath(@__DIR__, "..", "Test_Problem_Matrices")

A_system = readdlm(
    joinpath(matrixdir, "matrixN20x2.csv"),
    ',',
    Float64,
)

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

# Spectrum
S = ComplexF64.(eigvals(A_system))
dt_inner = 5 * 1e-4


# ============================================================
# Table 4.51: Maximal Stable Step Sizes of PFE and PRK4
# ============================================================
println()
println("Table 4.51: Maximal stable step sizes of PFE and PRK4.")

for K in 7:9
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
        dt_max_search = 0.5,
    )
end



# ============================================================
# Table 4.52: Largest successfully constructed stage numbers
# for the generalized two-scale spectrum
# ============================================================

# Chebyshev basis with 3S* formulation
matdir_3S_chebyshev  = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_3Sstar", "Test_3Sstar_chebychev")
algs_3S_chebyshev = load_algs_auto(matdir_3S_chebyshev; method_type = "3Sstar");
dt_map_3S_chebyshev = load_dt_map_from_polyopt_data(matdir_3S_chebyshev)

println()
println("Table 4.52: Chebyshev basis with 3S* formulation.")
data_stab = compare_stability_on_spectrum(algs_3S_chebyshev;
    S = S,
    dt_map = dt_map_3S_chebyshev,
    stab_tol = 1e-8
);

# Chebyshev basis with 3S* formulation and monotone c
matdir_3S_chebyshev_monotone  = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_3Sstar", "Test_3Sstar_chebychev,monotone")
algs_3S_chebyshev_monotone = load_algs_auto(matdir_3S_chebyshev_monotone; method_type = "3Sstar");
dt_map_3S_chebyshev_monotone = load_dt_map_from_polyopt_data(matdir_3S_chebyshev_monotone)

println()
println("Table 4.52: Chebyshev basis with 3S* formulation and monotone c.")
data_stab = compare_stability_on_spectrum(algs_3S_chebyshev_monotone;
    S = S,
    dt_map = dt_map_3S_chebyshev_monotone,
    stab_tol = 1e-8
);




# Chebyshev basis with 2N formulation
matdir_2N_chebyshev = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_2N", "Test_2N_chebychev")
algs_2N_chebyshev = load_algs_auto(matdir_2N_chebyshev; method_type = "2N");
dt_map_2N_chebyshev = load_dt_map_from_polyopt_data(matdir_2N_chebyshev)

println()
println("Table 4.52: Chebyshev basis with 2N formulation.")
data_stab = compare_stability_on_spectrum(algs_2N_chebyshev;
    S = S,
    dt_map = dt_map_2N_chebyshev,
    stab_tol = 1e-8
);

# Chebyshev basis with 2N formulation and monotone c
matdir_2N_monotone = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_2N", "Test_2N_chebychev_monotone")
algs_2N_monotone = load_algs_auto(matdir_2N_monotone; method_type = "2N");
dt_map_2N_monotone = load_dt_map_from_polyopt_data(matdir_2N_monotone)

println()
println("Table 4.52: Chebyshev basis with 2N formulation and monotone c.")
data_stab = compare_stability_on_spectrum(algs_2N_monotone;
    S = S,
    dt_map = dt_map_2N_monotone,
    stab_tol = 1e-8
);

# Chebyshev interpolation with 3S* formulation and interpolation conditions located near the fast cluster
matdir_3Sstar_fast_focused = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_3Sstar", "Test_3Sstar_chebyinterpol","3Sstar_fast_focused")
algs_3Sstar_fast_focused = load_algs_auto(matdir_3Sstar_fast_focused; method_type = "3Sstar");
dt_map_3Sstar_fast_focused = load_dt_map_from_polyopt_data(matdir_3Sstar_fast_focused)

println()
println("Table 4.52: Chebyshev interpolation with 3S* formulation and interpolation conditions located near the fast cluster.")
data_stab = compare_stability_on_spectrum(algs_3Sstar_fast_focused;
    S = S,
    dt_map = dt_map_3Sstar_fast_focused,
    stab_tol = 1e-8
);

# Chebyshev interpolation with 3S* formulation and interpolation conditions spread over the full interval
matdir_3Sstar_full_interval = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_full_interval")
algs_3Sstar_full_interval = load_algs_auto(matdir_3Sstar_full_interval; method_type = "3Sstar");
dt_map_3Sstar_full_interval = load_dt_map_from_polyopt_data(matdir_3Sstar_full_interval)

println()
println("Table 4.52: Chebyshev interpolation with 3S* formulation and interpolation conditions spread over the full interval.")
data_stab = compare_stability_on_spectrum(algs_3Sstar_full_interval;
    S = S,
    dt_map = dt_map_3Sstar_full_interval,
    stab_tol = 1e-8
);



# Chebyshev interpolation with 3S* formulation and interpolation conditions divided between the slow and fast clusters
matdir_3Sstar_cluster_split = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_cluster_split")
algs_3Sstar_cluster_split = load_algs_auto(matdir_3Sstar_cluster_split; method_type = "3Sstar");
dt_map_3Sstar_cluster_split = load_dt_map_from_polyopt_data(matdir_3Sstar_cluster_split)

println()
println("Table 4.52: Chebyshev interpolation with 3S* formulation and interpolation conditions divided between the slow and fast clusters.")
data_stab = compare_stability_on_spectrum(algs_3Sstar_cluster_split;
    S = S,
    dt_map = dt_map_3Sstar_cluster_split,
    stab_tol = 1e-8
);


# ============================================================
# Table 4.53: Derivative-Based PolyOpt Order Conditions for
# Chebyshev Interpolation with Conditions Spread over the
# Full Interval
# ============================================================

matdir_3Sstar_full_interval_derivative = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_full_interval_derivative")
algs_3Sstar_full_interval_derivative = load_algs_auto(matdir_3Sstar_full_interval_derivative; method_type = "3Sstar");
dt_map_3Sstar_full_interval_derivative = load_dt_map_from_polyopt_data(matdir_3Sstar_full_interval_derivative)

println()
println(
    "Table 4.53: Derivative-based PolyOpt order conditions with interpolation conditions spread over the full interval.",
)
data_stab = compare_stability_on_spectrum(algs_3Sstar_full_interval_derivative;
    S = S,
    dt_map = dt_map_3Sstar_full_interval_derivative,
    stab_tol = 1e-8
);

# ============================================================
# Table 4.54: Optimized 3S* Methods with chebyshev Basis at Polyopt stepsize
# ============================================================

# Chebyshev basis with 3S* formulation
matdir_3S_chebyshev  = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-20_lamF=-2000", "Test_3Sstar", "Test_3Sstar_chebychev")
algs_3S_chebyshev = load_algs_auto(matdir_3S_chebyshev; method_type = "3Sstar");
dt_map_3S_chebyshev = load_dt_map_from_polyopt_data(matdir_3S_chebyshev)

println()
println("Table 4.54: Optimized 3S* Methods with chebyshev Basis at Polyopt stepsize.")
data_stab = compare_stability_on_spectrum(algs_3S_chebyshev;
    S = S,
    dt_map = dt_map_3S_chebyshev,
    stab_tol = 1e-8
);


# ============================================================
# Table 4.55: Comparison of Butcher and Low-Storage
# Stability Evaluations for the Selected Methods
# ============================================================

println()
println(
    "Table 4.55: Comparison of Butcher and low-storage stability evaluations for the selected methods.",
)

selected_methods_generalized = [
    "3Sstar-1-13_2026-06-30T16-53-41.txt",
    "3Sstar-2-15_2026-06-29T00-30-14.txt",
    "3Sstar-3-18_2026-06-28T04-25-30.txt",
    "3Sstar-4-20_2026-07-05T09-03-02.txt",
]

data_stab_generalized = print_stability_evaluation_table(
    matdir_3S_chebyshev,
    selected_methods_generalized,
    S;
    method_type = "3Sstar",
)


# ============================================================
# Selected Methods for Tables 4.56-4.57
# ============================================================

algs_selected, dt_map_selected = select_methods(
    algs_3S_chebyshev,
    dt_map_3S_chebyshev,
    selected_methods_generalized,
)

# Compute the maximal stable outer step sizes by bisection
K = 8
dt_max_PFE, _ = find_dt_max_cluster_PI(
    K = K,
    dt_inner = dt_inner,
    S = S,
    c_out = cout_PFE,
    A_out = Aout_PFE,
    b_out = bout_PFE,
    dt_max_search = 0.5,
    stab_tol = 0
);

dt_max_PRK4, _ = find_dt_max_cluster_PI(
    K = K,
    dt_inner = dt_inner,
    S = S,
    c_out = cout_RK4,
    A_out = Aout_RK4,
    b_out = bout_RK4,
    dt_max_search = 0.5,
    stab_tol = 0
);


# ============================================================
# Table 4.56: Error Comparison for T = 0.1067
# ============================================================
u0 = vcat(collect(1.0:20), zeros(20))
T = 0.1067

println()
println("Table 4.56, RK-Opt methods: Errors for T = 0.1067.")
data_ode_1 = compare_generalized_two_scale_ode(
    algs_selected;
    A_sys = A_system,
    S = S,
    u0 = u0,
    dt_map = dt_map_selected,
    tspan = (0.0, T),
    stab_tol = 1e-8,
)

println()
println("Table 4.56, PI methods: Errors for T = 0.1067.")
ode_data = compare_PI_methods_generalized_two_scale_ode(
    outer_tableaus = [
        (name = "PFE", c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 8, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 8, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    A_sys = A_system,
    S = S,
    u0 = u0,
    tspan = (0.0, T),
)






# ============================================================
# Table 4.57: Error Comparison for T = 50
# ============================================================
T = 50 

println()
println("Table 4.57, RK-Opt methods: Errors for T = 50.")
data_ode_1 = compare_generalized_two_scale_ode(
    algs_selected;
    A_sys = A_system,
    S = S,
    u0 = u0,
    dt_map = dt_map_selected,
    tspan = (0.0, T),
    stab_tol = 1e-8,
)

println()
println("Table 4.57, PI methods: Errors for T = 50.")
ode_data = compare_PI_methods_generalized_two_scale_ode(
    outer_tableaus = [
        (name = "PFE", c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 8, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 8, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    A_sys = A_system,
    S = S,
    u0 = u0,
    tspan = (0.0, T),
)



# ============================================================
# Table 4.58: Error of PFE and PRK4 with reduced Step Sizes
# ============================================================
dt_PFE  = 0.045
dt_PRK4 = 0.0872
T = 50

println()
println("Table 4.58: PI methods with reduced step sizes.")
ode_data = compare_PI_methods_generalized_two_scale_ode(
    outer_tableaus = [
        (name = "PFE", c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 8, dt_inner = dt_inner, dt = dt_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 8, dt_inner = dt_inner, dt = dt_PRK4, p = 4),
    ],
    A_sys = A_system,
    S = S,
    u0 = u0,
    tspan = (0.0, T),
)

