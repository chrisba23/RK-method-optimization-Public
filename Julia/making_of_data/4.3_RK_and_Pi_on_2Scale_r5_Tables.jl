# Reproduces the data for the two-scale problem with r=5  in Chapter 4.3
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

# Two-scale cluster spectrum with r_fast = 5
S = make_eig_clusters(
    α = α,
    ε = ε,
    kappa_slow = 0.005,
    kappa_fast = 0.005
)

# ============================================================
# Table 4.12: Maximal stable step sizes of PFE and PRK4
# ============================================================
println()
println("Table 4.12")
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
        dt_max_search = 5.0
    )
end




# ============================================================
# Table 4.18: Optimized 3S* Methods
# ============================================================
method_type = "3Sstar"

# Methods using step sizes from the additional grid search
matdir_grid = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=5_linsearch", "3Sstar_fast_focused")
algs_grid  = load_algs_auto(matdir_grid; method_type = method_type)
dt_map_grid  = load_dt_map_from_polyopt_data(matdir_grid)

println()
println("Table 4.18: Optimized 3S* methods using the additional grid-search step sizes.")
data_stab = compare_stability_on_spectrum(algs_grid;
    S = S,
    dt_map = dt_map_grid,
    stab_tol = 1e-8
);

# Methods using the original PolyOpt step sizes
matdir_polyopt  = joinpath(@__DIR__, "..", "..", "Matlab", "two_scale_lamS=-1_lamF=-1000", "Test_3Sstar", "Test_3Sstar_chebyinterpol", "2-Scale,r=5", "3Sstar_fast_focused")
algs_polyopt  = load_algs_auto(matdir_polyopt; method_type = method_type)
dt_map_polyopt  = load_dt_map_from_polyopt_data(matdir_polyopt)

println()
println("Table 4.18: Optimized 3S* methods at the original PolyOpt step sizes.")
data_stab = compare_stability_on_spectrum(algs_polyopt;
    S = S,
    dt_map = dt_map_polyopt,
    stab_tol = 1e-8
);




# ============================================================
# Selected Methods for Tables 4.19-4.21
# ============================================================
# Combine methods from both folders for the error comparisons
algs_all = vcat(algs_grid, algs_polyopt)
dt_map_all = merge(dt_map_grid, dt_map_polyopt)

selected_files = [
    "3Sstar-1-6_2026-04-04T22-26-56.txt",
    "3Sstar-3-8_2026-05-09T14-43-33.txt",
    "3Sstar-4-9_2026-08-01T13-57-04.txt",
    "3Sstar-2-11_2026-07-22T14-14-33.txt",
]

algs_selected, dt_map_selected = select_methods(
    algs_all,
    dt_map_all,
    selected_files,
)

# Compute the maximal stable outer step sizes by bisection
K = 2
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
# Table 4.19: Error Comparison for T = 27.4
# ============================================================

T = 27.4
println()
println("Table 4.19, RK-Opt methods: Errors for T = 27.4.")
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
println("Table 4.19, PI methods: Errors for T = 27.4.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 2, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 2, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    problem = :cluster,
    α = α,
    ε = ε,
    S = S,
    tspan = (0.0, T)
)
# ============================================================
# Table 4.20: Error Comparison for T = 202
# ============================================================

T = 202
println()
println("Table 4.20/Table 4.21, RK-Opt methods: Errors for T = 202.")
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
println("Table 4.20, PI methods: Errors for T = 202.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 2, dt_inner = dt_inner, dt = dt_max_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 2, dt_inner = dt_inner, dt = dt_max_PRK4, p = 4),
    ],
    problem = :cluster,
    α = α,
    ε = ε,
    S = S,
    tspan = (0.0, T)
)

# ============================================================
# Table 4.21: Component-wise Errors and Stability Values
# ============================================================

# The component-wise errors are computed in Table 4.20
println()
println("Table 4.21: Stability values on the slow and fast clusters.")

rk_cluster_stab = clusterwise_stability_RK(
    algs_selected;
    S = S,
    dt_map = dt_map_selected
)



# ============================================================
# Table 4.22: PI Methods with Reduced Step Sizes
# ============================================================
dt_PFE  = 1.2
dt_PRK4 = 1.11
T = 202
println()
println("Table 4.22: Errors of PFE and PRK4 with reduced outer step sizes.")
ode_data = compare_PI_methods_ode(
    outer_tableaus = [
        (name = "PFE",  c = cout_PFE, A = Aout_PFE, b = bout_PFE,
         K = 2, dt_inner = dt_inner, dt = dt_PFE, p = 1),

        (name = "PRK4", c = cout_RK4, A = Aout_RK4, b = bout_RK4,
         K = 2, dt_inner = dt_inner, dt = dt_PRK4, p = 4),
    ],
    problem = :cluster,
    α = α,
    ε = ε,
    S = S,
    tspan = (0.0, T)
)

