# Reproduces the data for the hyperbolic problem in Chapter 4.5
include(joinpath(@__DIR__, "..", "RK_Data_Functions.jl"))
include(joinpath(@__DIR__, "..", "Pi_Methods_Data.jl"))


# ============================================================
# Load Problem Matrices
# ============================================================

matrixdir = joinpath(@__DIR__, "..", "Test_Problem_Matrices")

A_eps1e3 = readdlm(
    joinpath(matrixdir, "matrix_eps10e-3.csv"),
    ',',
    Float64,
)

A_eps1e4 = readdlm(
    joinpath(matrixdir, "matrix_eps10e-4.csv"),
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



# ============================================================
# Hyperbolic problem with ε = 1e-3
# ============================================================
S1 = ComplexF64.(eigvals(A_eps1e3))
λfast_eps1e3 = -1299.928
dt_inner_eps1e3 = -1 / λfast_eps1e3

# ============================================================
# Table 4.41: Maximal Stable Step Sizes of PFE and PRK4
# for ε = 1e-3
# ============================================================

println()
println("Table 4.41: Maximal stable step sizes of PFE and PRK4 for ε = 1e-3.")
for K in 1:3
    println()
    println("K = $K.")


    outer_tableaus = [
        (
            name = "PFE",
            c = cout_PFE,
            A = Aout_PFE,
            b = bout_PFE,
            K = K,
            dt_inner = dt_inner_eps1e3,
        ),
        (
            name = "PRK4",
            c = cout_RK4,
            A = Aout_RK4,
            b = bout_RK4,
            K = K,
            dt_inner = dt_inner_eps1e3,
        ),
    ]

    compare_PI_methods_on_cluster(
        outer_tableaus = outer_tableaus,
        S = S1,
        dt_max_search = 0.05,
    )
end




# ============================================================
# Table 4.42: Largest successfully constructed stage numbers
# for the moment-model spectrum with ε = 1e-3
# ============================================================
# Chebyshev basis with 3S* formulation
matdir_eps1e3_3S_chebyshev  = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-1299.93_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyshev")
algs_eps1e3_3S_chebyshev = load_algs_auto(matdir_eps1e3_3S_chebyshev; method_type = "3Sstar");
dt_map_eps1e3_3S_chebyshev = load_dt_map_from_polyopt_data(matdir_eps1e3_3S_chebyshev)

println()
println("Table 4.42: Chebyshev basis with 3S* formulation.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_3S_chebyshev;
    S = S1,
    dt_map = dt_map_eps1e3_3S_chebyshev,
    stab_tol = 1e-8
);


# Chebyshev basis with 3S* formulation and monotone c
matdir_eps1e3_3S_chebyshev_monotone  = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-1299.93_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyshev,monotone")
algs_eps1e3_3S_chebyshev_monotone = load_algs_auto(matdir_eps1e3_3S_chebyshev_monotone; method_type = "3Sstar");
dt_map_eps1e3_3S_chebyshev_monotone = load_dt_map_from_polyopt_data(matdir_eps1e3_3S_chebyshev_monotone)

println()
println("Table 4.42: Chebyshev basis with 3S* formulation and monotone c.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_3S_chebyshev_monotone;
    S = S1,
    dt_map = dt_map_eps1e3_3S_chebyshev_monotone,
    stab_tol = 1e-8
);



# Chebyshev basis with 2N formulation
matdir_eps1e3_2N_chebyshev = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-1299.93_R=299.928", "Test_2N", "Test_2N_chebyshev")
algs_eps1e3_2N_chebyshev = load_algs_auto(matdir_eps1e3_2N_chebyshev; method_type = "2N");
dt_map_eps1e3_2N_chebyshev = load_dt_map_from_polyopt_data(matdir_eps1e3_2N_chebyshev)

println()
println("Table 4.42: Chebyshev basis with 2N formulation.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_2N_chebyshev;
    S = S1,
    dt_map = dt_map_eps1e3_2N_chebyshev,
    stab_tol = 1e-8
);

# Chebyshev basis with 2N formulation and nondecreasing stage nodes
matdir_eps1e3_2N_monotone = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-1299.93_R=299.928", "Test_2N", "Test_2N_chebyshev,monotone")
algs_eps1e3_2N_monotone = load_algs_auto(matdir_eps1e3_2N_monotone; method_type = "2N");
dt_map_eps1e3_2N_monotone = load_dt_map_from_polyopt_data(matdir_eps1e3_2N_monotone)

println()
println("Table 4.42: Chebyshev basis with 2N formulation and nondecreasing stage nodes.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_2N_monotone;
    S = S1,
    dt_map = dt_map_eps1e3_2N_monotone,
    stab_tol = 1e-8
);


# Chebyshev interpolation with 3S* formulation and interpolation conditions located near the fast cluster
matdir_eps1e3_interp_fast  = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-1299.93_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_fast_focused")
algs_eps1e3_interp_fast = load_algs_auto(matdir_eps1e3_interp_fast; method_type = "3Sstar");
dt_map_eps1e3_interp_fast = load_dt_map_from_polyopt_data(matdir_eps1e3_interp_fast)

println()
println("Table 4.42: Chebyshev interpolation with 3S* formulation and interpolation conditions located near the fast cluster.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_interp_fast;
    S = S1,
    dt_map = dt_map_eps1e3_interp_fast,
    stab_tol = 1e-8
);


# Chebyshev interpolation with 3S* formulation and interpolation conditions spread over the full interval
matdir_eps1e3_interp_full = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-1299.93_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_full_interval")
algs_eps1e3_interp_full = load_algs_auto(matdir_eps1e3_interp_full; method_type = "3Sstar");
dt_map_eps1e3_interp_full = load_dt_map_from_polyopt_data(matdir_eps1e3_interp_full)

println()
println("Table 4.42: Chebyshev interpolation with 3S* formulation and interpolation conditions spread over the full interval.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_interp_full;
    S = S1,
    dt_map = dt_map_eps1e3_interp_full,
    stab_tol = 1e-8
);

# Chebyshev interpolation with 3S* formulation and interpolation conditions divided between the slow and fast clusters
matdir_eps1e3_interp_cluster = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-1299.93_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_cluster_split")
algs_eps1e3_interp_cluster = load_algs_auto(matdir_eps1e3_interp_cluster; method_type = "3Sstar");
dt_map_eps1e3_interp_cluster = load_dt_map_from_polyopt_data(matdir_eps1e3_interp_cluster)

println()
println("Table 4.42: Chebyshev interpolation with 3S* formulation and interpolation conditions divided between the slow and fast clusters.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_interp_cluster;
    S = S1,
    dt_map = dt_map_eps1e3_interp_cluster,
    stab_tol = 1e-8
);




# ============================================================
# Table 4.43: Derivative-Based PolyOpt Order Conditions for
# Chebyshev Interpolation with Cluster-Split Matching Points
# for ε = 1e-3
# ============================================================

matdir_eps1e3_cluster_split_derivative = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-1299.93_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_cluster_split_derivative")
algs_eps1e3_cluster_split_derivative = load_algs_auto(matdir_eps1e3_cluster_split_derivative; method_type = "3Sstar");
dt_map_eps1e3_cluster_split_derivative = load_dt_map_from_polyopt_data(matdir_eps1e3_cluster_split_derivative)

println()
println(
    "Table 4.43: Derivative-based PolyOpt order conditions with interpolation conditions divided between the slow and fast clusters.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_cluster_split_derivative;
    S = S1,
    dt_map = dt_map_eps1e3_cluster_split_derivative,
    stab_tol = 1e-8
);




# ============================================================
# Table 4.44: Optimized 3S* methods with Chebyshev basis
# ============================================================

println()
println("Table 4.44: Optimized 3S* methods with Chebyshev basis.")
data_stab = compare_stability_on_spectrum(algs_eps1e3_3S_chebyshev;
    S = S1,
    dt_map = dt_map_eps1e3_3S_chebyshev,
    stab_tol = 1e-8
);



# ============================================================
# Values below Table 4.44: Comparison of Butcher and
# Low-Storage Stability Evaluations for the Selected Methods
# for ε = 1e-3
# ============================================================

println()
println("Values below Table 4.44: Comparison of Butcher and low-storage stability evaluations for the selected methods with ε = 1e-3.")

selected_methods_eps1e3 = [
    "3Sstar-1-10_2026-07-02T09-42-57.txt",
    "3Sstar-2-20_2026-06-30T15-34-45.txt",
    "3Sstar-3-19_2026-06-26T02-29-01.txt",
    "3Sstar-4-21_2026-07-04T19-40-26.txt",
]

data_stab_eps1e3 = print_stability_evaluation_table(
    matdir_eps1e3_3S_chebyshev,
    selected_methods_eps1e3,
    S1;
    method_type = "3Sstar",
)


# ============================================================
# Table 4.45: Selected Methods for Comparison with PI Methods
# ============================================================

algs_selected_eps1e3, dt_map_selected_eps1e3 = select_methods(
    algs_eps1e3_3S_chebyshev,
    dt_map_eps1e3_3S_chebyshev,
    selected_methods_eps1e3,
)

println()
println(
    "Table 4.45: Selected 3S* methods for comparison with PI methods for ε = 1e-3.",
)

data_stab_selected_eps1e3 = compare_stability_on_spectrum(
    algs_selected_eps1e3;
    S = S1,
    dt_map = dt_map_selected_eps1e3,
    stab_tol = 1e-8,
)




# ============================================================
# Hyperbolic problem with ε = 1e-4
# ============================================================
S2 = ComplexF64.(eigvals(A_eps1e4))
λfast_eps1e4 = -10299.928
dt_inner_eps1e4 = -1 / λfast_eps1e4


# ============================================================
# Table 4.46: Maximal Stable Step Sizes of PFE and PRK4
# for ε = 1e-4
# ============================================================


println()
println("Table 4.46: Maximal stable step sizes of PFE and PRK4 for ε = 1e-4.")
for K in 1:2
    println("\nK = $K")

    outer_tableaus = [
        (
            name = "PFE",
            c = cout_PFE,
            A = Aout_PFE,
            b = bout_PFE,
            K = K,
            dt_inner = dt_inner_eps1e4,
        ),
        (
            name = "PRK4",
            c = cout_RK4,
            A = Aout_RK4,
            b = bout_RK4,
            K = K,
            dt_inner = dt_inner_eps1e4,
        ),
    ]

    compare_PI_methods_on_cluster(
        outer_tableaus = outer_tableaus,
        S = S2,
        dt_max_search = 0.05,
    )
end







# ============================================================
# Table 4.47: Largest successfully constructed stage numbers
# for the moment-model spectrum with ε = 1e-4
# ============================================================

# Chebyshev basis with 3S* formulation
matdir_eps1e4_3S_chebyshev  = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyshev")
method_type = "3Sstar"
algs_eps1e4_3S_chebyshev = load_algs_auto(matdir_eps1e4_3S_chebyshev; method_type = "3Sstar");
dt_map_eps1e4_3S_chebyshev = load_dt_map_from_polyopt_data(matdir_eps1e4_3S_chebyshev)

println()
println("Table 4.47: Chebyshev basis with 3S* formulation.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_3S_chebyshev;
    S = S2,
    dt_map = dt_map_eps1e4_3S_chebyshev,
    stab_tol = 1e-8
);


# Chebyshev basis with 3S* formulation and monotone c
matdir_eps1e4_3S_chebyshev_monotone  = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyshev,monotone")
method_type = "3Sstar"
algs_eps1e4_3S_chebyshev_monotone = load_algs_auto(matdir_eps1e4_3S_chebyshev_monotone; method_type = "3Sstar");
dt_map_eps1e4_3S_chebyshev_monotone = load_dt_map_from_polyopt_data(matdir_eps1e4_3S_chebyshev_monotone)

println()
println("Table 4.47: Chebyshev basis with 3S* formulation and monotone c.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_3S_chebyshev_monotone;
    S = S2,
    dt_map = dt_map_eps1e4_3S_chebyshev_monotone,
    stab_tol = 1e-8
);



# Chebyshev basis with 2N formulation
matdir_eps1e4_2N_chebyshev = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_2N", "Test_2N_chebyshev")
algs_eps1e4_2N_chebyshev = load_algs_auto(matdir_eps1e4_2N_chebyshev; method_type = "2N");
dt_map_eps1e4_2N_chebyshev = load_dt_map_from_polyopt_data(matdir_eps1e4_2N_chebyshev)

println()
println("Table 4.47: Chebyshev basis with 2N formulation.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_2N_chebyshev;
    S = S2,
    dt_map = dt_map_eps1e4_2N_chebyshev,
    stab_tol = 1e-8
);


# Chebyshev basis with 2N formulation and monotone c
matdir_eps1e4_2N_monotone = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_2N", "Test_2N_chebyshev,monotone")
algs_eps1e4_2N_monotone = load_algs_auto(matdir_eps1e4_2N_monotone; method_type = "2N");
dt_map_eps1e4_2N_monotone = load_dt_map_from_polyopt_data(matdir_eps1e4_2N_monotone)

println()
println("Table 4.47: Chebyshev basis with 2N formulation and monotone c.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_2N_monotone;
    S = S2,
    dt_map = dt_map_eps1e4_2N_monotone,
    stab_tol = 1e-8
);



# Chebyshev interpolation with 3S* formulation and interpolation conditions located near the fast cluster
matdir_eps1e4_interp_fast  = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_fast_focused")
algs_eps1e4_interp_fast = load_algs_auto(matdir_eps1e4_interp_fast; method_type = "3Sstar");
dt_map_eps1e4_interp_fast = load_dt_map_from_polyopt_data(matdir_eps1e4_interp_fast)

println()
println("Table 4.47: Chebyshev interpolation with 3S* formulation and interpolation conditions located near the fast cluster.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_interp_fast;
    S = S2,
    dt_map = dt_map_eps1e4_interp_fast,
    stab_tol = 1e-8
);

# Chebyshev interpolation with 3S* formulation and interpolation conditions spread over the full interval
matdir_eps1e4_interp_full  = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_full_interval")
algs_eps1e4_interp_full = load_algs_auto(matdir_eps1e4_interp_full; method_type = "3Sstar");
dt_map_eps1e4_interp_full = load_dt_map_from_polyopt_data(matdir_eps1e4_interp_full)

println()
println("Table 4.47: Chebyshev interpolation with 3S* formulation and interpolation conditions spread over the full interval.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_interp_full;
    S = S2,
    dt_map = dt_map_eps1e4_interp_full,
    stab_tol = 1e-8
);


# Chebyshev interpolation with 3S* formulation and interpolation conditions divided between the slow and fast clusters
matdir_eps1e4_cluster_split = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_cluster_split")
algs_eps1e4_cluster_split = load_algs_auto(matdir_eps1e4_cluster_split; method_type = "3Sstar");
dt_map_eps1e4_cluster_split = load_dt_map_from_polyopt_data(matdir_eps1e4_cluster_split)

println()
println("Table 4.47: Chebyshev interpolation with 3S* formulation and interpolation conditions divided between the slow and fast clusters.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_cluster_split;
    S = S2,
    dt_map = dt_map_eps1e4_cluster_split,
    stab_tol = 1e-8
);




# ============================================================
# Table 4.48: Derivative-Based PolyOpt Order Conditions for
# Chebyshev Interpolation with Conditions Divided Between the
# Slow and Fast Clusters for ε = 1e-4
# ============================================================

matdir_eps1e4_cluster_split_derivative = joinpath(@__DIR__, "..", "..", "Matlab", "hyp_upwind_model", "lamS=-299.928_lamF=-10299.9_R=299.928", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "3Sstar_cluster_split_derivative")
algs_eps1e4_cluster_split_derivative = load_algs_auto(matdir_eps1e4_cluster_split_derivative; method_type = "3Sstar");
dt_map_eps1e4_cluster_split_derivative = load_dt_map_from_polyopt_data(matdir_eps1e4_cluster_split_derivative)

println()
println(
    "Table 4.48: Derivative-based PolyOpt order conditions with interpolation conditions divided between the slow and fast clusters.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_cluster_split_derivative;
    S = S2,
    dt_map = dt_map_eps1e4_cluster_split_derivative,
    stab_tol = 1e-8
);


# ============================================================
# Table 4.49: Optimized 3S* methods with Chebyshev basis
# ============================================================
println()
println("Table 4.49: Optimized 3S* methods with Chebyshev basis.")
data_stab = compare_stability_on_spectrum(algs_eps1e4_3S_chebyshev;
    S = S2,
    dt_map = dt_map_eps1e4_3S_chebyshev,
    stab_tol = 1e-8
);




# ============================================================
# Values below Table 4.49: Comparison of Butcher and
# Low-Storage Stability Evaluations for the Selected Methods
# for ε = 1e-4
# ============================================================

println()
println(
    "Values below Table 4.49: Comparison of Butcher and low-storage stability evaluations for the selected methods with ε = 1e-4.",
)

selected_methods_eps1e4 = [
    "3Sstar-1-10_2026-06-28T23-27-09.txt",
    "3Sstar-2-12_2026-06-30T13-58-54.txt",
    "3Sstar-3-15_2026-05-13T00-07-09.txt",
    "3Sstar-4-15_2026-05-12T16-29-07.txt",
]

data_stab_eps1e4 = print_stability_evaluation_table(
    matdir_eps1e4_3S_chebyshev,
    selected_methods_eps1e4,
    S2;
    method_type = "3Sstar",
)


# ============================================================
# Table 4.50: Selected Methods for Comparison with PI Methods
# ============================================================

algs_selected_eps1e4, dt_map_selected_eps1e4 = select_methods(
    algs_eps1e4_3S_chebyshev,
    dt_map_eps1e4_3S_chebyshev,
    selected_methods_eps1e4,
)

println()
println(
    "Table 4.50: Selected 3S* methods for comparison with PI methods for ε = 1e-4.",
)

data_stab_selected_eps1e4 = compare_stability_on_spectrum(
    algs_selected_eps1e4;
    S = S2,
    dt_map = dt_map_selected_eps1e4,
    stab_tol = 1e-8,
)