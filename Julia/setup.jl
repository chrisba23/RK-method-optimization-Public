# ============================================================
# Shared helper functions and setups
# ============================================================

import Pkg
Pkg.activate(@__DIR__)
#Pkg.instantiate() run once

using DiffEqBase
using OrdinaryDiffEq
using MAT
using Plots
using PrettyTables: pretty_table, ft_printf
using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles

include(joinpath(@__DIR__, "..", "Extern-Repos", "2025_hyp_RK_dev-main", "code_julia", "code.jl"))

# Create the slow and fast eigenvalue clusters for the two-scale test problem.
function make_eig_clusters(; α = 1.0, ε = 1e-3,
    kappa_slow = 0.0, kappa_fast = 0.0,
    M_slow = 300, M_fast = 300)

    λs = -α
    λf = -1 / ε

    if kappa_slow == 0.0
        Dslow = [complex(λs)]
    else
        r1 = kappa_slow * abs(λs)
        θ1 = range(0, 2π, length=M_slow)
        Dslow = λs .+ r1*cos.(θ1) .+ im*r1*sin.(θ1)
    end

    if kappa_fast == 0.0
        Dfast = [complex(λf)]
    else
        r2 = kappa_fast * abs(λf)
        θ2 = range(0, 2π, length=M_fast)
        Dfast = λf .+ r2*cos.(θ2) .+ im*r2*sin.(θ2)
    end

    return vcat(Dslow, Dfast)
end


# Evaluate the stability function from a Butcher tableau.
function R_general(A, b, z)
    s = length(b)
    e = ones(s)

    stage_values = (I - z * A) \ e

    return 1 + z * sum(b .* stage_values)
end

# Evaluate the stability function on all entries of Z.
function R_grid(A, b, Z)
    vals = [R_general(A, b, z) for z in vec(Z)]
    return reshape(vals, size(Z))
end



# Exact solution for the two-scale ODE test problem
function exact_two_scale(t; α = 1.0, ε = 1e-3, u10 = 1.0, u20 = 0.0)
    u1 = u10 * exp(-α * t)
    b  = 1 - α * ε
    u2 = (u10 / b) * exp(-α * t) + (u20 - u10 / b) * exp(-t / ε)
    return u1, u2
end

# Solve the two-scale ODE test problem
function solve_two_scale(alg; dt,
    α = 1.0,
    ε = 1e-3,
    tspan = (0.0, 1.0),
    u10 = 1.0,
    u20 = 0.0)

    function f!(du, u, p, t)
        du[1] = -α * u[1]
        du[2] = (u[1] - u[2]) / ε
    end

    prob = ODEProblem(f!, [u10, u20], tspan)
    return solve(prob, alg; dt = dt, save_everystep = false, adaptive = false)
end


# Exact solution for the two-scale Cluster ODE test problem
function exact_cluster(t, S; u0 = ones(ComplexF64, length(S)))
    return exp.(S .* t) .* u0
end

# Solve the two-scale Cluster ODE test problem
function solve_cluster_ode(alg; S, dt,
    tspan = (0.0, 1.0),
    u0 = ones(ComplexF64, length(S)))

    function f!(du, u, p, t)
        du .= S .* u
    end

    prob = ODEProblem(f!, u0, tspan)

    return solve(prob, alg; dt = dt, save_everystep = false, adaptive = false)
end

# Solve the generalized two-scale ODE test problem
function solve_generalized_two_scale(alg; A_sys, dt,
    tspan = (0.0, 1.0),
    u0 = ones(ComplexF64, size(A_sys, 1)))

    function f!(du, u, p, t)
        mul!(du, A_sys, u)
    end

    prob = ODEProblem(f!, u0, tspan)

    return solve(prob, alg;
        dt = dt,
        save_everystep = false,
        adaptive = false
    )
end


# Construct the Butcher tableau of a Projective Runge-Kutta method.
function create_PRK_tableau(c, A, b, K, λ)
    n_outer_stages = length(c)
    n_inner_stages = K + 1
    n_total_stages = n_inner_stages * n_outer_stages

    c_prk = zeros(n_total_stages)
    A_prk = zeros(n_total_stages, n_total_stages)
    b_prk = zeros(n_total_stages)

    e = ones(n_inner_stages)

    e_Kplus1 = zeros(n_inner_stages)
    e_Kplus1[end] = 1.0

    for s in 1:n_outer_stages
        idx_s = ((s - 1) * n_inner_stages + 1):(s * n_inner_stages)

        # Stage nodes:
        # λ(0, ..., K) for the first outer stage and
        # c_s + λ(0, ..., K) for all later outer stages.
        c_base = (s == 1) ? 0.0 : c[s]
        c_prk[idx_s] = c_base .+ λ .* (0:K)

        # Diagonal Forward Euler block A_λ for the inner stages.
        A_λ =
            λ * LowerTriangular(ones(n_inner_stages, n_inner_stages)) -
            λ * I(n_inner_stages)

        A_prk[idx_s, idx_s] = A_λ

        # Coupling to all previous outer stages.
        if s > 1
            for l in 1:(s - 1)
                idx_l = ((l - 1) * n_inner_stages + 1):(l * n_inner_stages)

                a_tilde_sl =
                    ((c[s] / λ) - n_inner_stages) *
                    A[s, l] / c[s]

                # The first coupling block also contains the contribution
                # of the initial inner Forward Euler stages.
                if l == 1
                    Λ =
                        λ * (
                            e * e' +
                            a_tilde_sl * e * e_Kplus1'
                        )
                else
                    Λ = λ * a_tilde_sl * e * e_Kplus1'
                end

                A_prk[idx_s, idx_l] = Λ
            end
        end

        # Construct the weight block belonging to outer stage s.
        if s == 1
            b_prk[idx_s] =
                λ * e +
                (1 - n_inner_stages * λ) *
                b[1] *
                e_Kplus1
        else
            b_prk[idx_s] =
                (1 - n_inner_stages * λ) *
                b[s] *
                e_Kplus1
        end
    end

    return c_prk, A_prk, b_prk
end



# Split the spectrum into slow and fast clusters
function split_slow_fast(S)

    S_vec = ComplexF64.(vec(S))

    perm = sortperm(real.(S_vec); rev = true)
    S_sorted = S_vec[perm]

    split = argmax(abs.(diff(real.(S_sorted))))
    
    slow_idx = sort(perm[1:split])
    fast_idx = sort(perm[(split + 1):end])

    return S_vec[slow_idx], S_vec[fast_idx]
end



# Evaluate the stability function in 3S* low-storage form.
function R_3Sstar(beta, gamma1, gamma2, gamma3, delta, Z)
    m = length(delta)
    T = eltype(Z)

    S3 = ones(T, size(Z))
    S1 = ones(T, size(Z))
    S2 = zeros(T, size(Z))

    for i in 2:(m + 1)
        S2 .= S2 .+ delta[i - 1] .* S1

        S1 .= gamma1[i] .* S1 .+
             gamma2[i] .* S2 .+
             gamma3[i] .* S3 .+
             beta[i - 1] .* (Z .* S1)
    end

    return S1
end



# Extract low-storage coefficients from 3S* method
function extract_3Sstar_params(alg)
    beta_mat = alg.coefficients.β
    gamma1 = alg.coefficients.γ1
    gamma2 = alg.coefficients.γ2
    gamma3 = alg.coefficients.γ3
    delta = alg.coefficients.δ
    
    s = size(beta_mat, 1)
    beta = [beta_mat[i, i-1] for i in 2:s]
    
    gamma1 = vec(gamma1)
    gamma2 = vec(gamma2)
    gamma3 = vec(gamma3)
    delta = vec(delta)
    
    return beta, gamma1, gamma2, gamma3, delta
end

# Extract low-storage coefficients from 2N method
function extract_2N_params(alg)
    A2N = vec(alg.coefficients.A2N)
    B2N = vec(alg.coefficients.B2N)
    return A2N, B2N
end

# Evaluate the stability function in 2N low-storage form.
function R_2N(A2N, B2N, Z)
    T = eltype(Z)

    S1 = ones(T, size(Z))
    S2 = zeros(T, size(Z))

    for i in eachindex(B2N)
        S2 .= A2N[i] .* S2 .+ Z .* S1
        S1 .= S1 .+ B2N[i] .* S2
    end

    return S1
end

# Universal low-storage evaluation
function R_lowstorage(alg, method_type, Z)
    if method_type == "3Sstar"
        beta, g1, g2, g3, delta = extract_3Sstar_params(alg)
        return R_3Sstar(beta, g1, g2, g3, delta, Z)
    elseif method_type == "2N"
        A2N, B2N = extract_2N_params(alg)
        return R_2N(A2N, B2N, Z)
    else
        return R_grid(alg.coefficients.A, alg.coefficients.b, Z)
    end
end



# Load the saved MATLAB data of the corosponding RK method.
function load_experiment(matdir, method_file)
    base, _ = splitext(method_file)
    datafile = joinpath(matdir, "polyopt_data", base * "_plot_data.mat")

    isfile(datafile) || error("Plot data file not found: $datafile")

    vars = matread(datafile)
    return vars["experiment"]
end

# Load one RK method from its text file.
function load_method(matdir, method_file, method_type)
    rkfile = joinpath(matdir, method_file)
    isfile(rkfile) || error("RK method file not found: $rkfile")

    if method_type == "3Sstar"
        return RK3SFromFile(filename = rkfile)
    elseif method_type == "2N"
        return RK2NFromFile(filename = rkfile)
    elseif method_type == "ERK"
        return ExplicitRKFromFile(filename = rkfile)
    else
        error("Unknown method_type: $method_type")
    end
end

# Compute stability polynomial coefficients from a Butcher tableau.
function stability_poly_coeffs(A, b)
    s = length(b)
    T = promote_type(eltype(A), eltype(b))

    e = ones(T, s)
    coeffs = zeros(T, s + 1)
    coeffs[1] = one(T)

    v = copy(e)

    for j in 1:s
        coeffs[j + 1] = sum(b .* v)
        v = A * v
    end

    return coeffs
end


# Evaluate the Chebyshev interpolation polynomial in barycentric form.
function bary_eval_cheb(z_nodes, f_nodes, z_eval)
    z_nodes = vec(z_nodes)
    f_nodes = vec(f_nodes)
    z_eval  = vec(z_eval)

    N = length(z_nodes) - 1

    w = ones(ComplexF64, N + 1)
    w[1] *= 0.5
    w[end] *= 0.5

    for j in 0:N
        w[j+1] *= (-1)^j
    end

    p = Vector{ComplexF64}(undef, length(z_eval))

    for k in eachindex(z_eval)
        z = z_eval[k]
        idx = findfirst(abs.(z .- z_nodes) .< 1e-14)

        if idx !== nothing
            p[k] = f_nodes[idx]
        else
            tmp = w ./ (z .- z_nodes)
            p[k] = sum(tmp .* f_nodes) / sum(tmp)
        end
    end

    return p
end


# Evaluate the saved PolyOpt reference stability polynomial.
function eval_polyopt_reference(experiment, Z)
    basis = String(experiment["basis"])

    if basis == "chebyshev"
        a = vec(experiment["a"])
        return map(z -> evalpoly(z, a), Z)

    elseif basis == "chebinterp"
        z_nodes_full = vec(experiment["z_nodes_full"])
        R_vals_full  = vec(experiment["R_vals_full"])

        vals = bary_eval_cheb(z_nodes_full, R_vals_full, vec(Z))
        return reshape(vals, size(Z))

    else
        error("Unknown basis: $basis")
    end
end