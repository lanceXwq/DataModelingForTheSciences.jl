### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 1398e536-09ff-11ec-177f-6b6668e95c95
begin
    import Pkg
    Pkg.activate("")
    using Distributions, Plots, LaTeXStrings
end

# ╔═╡ 21524504-7968-4a8d-b99f-1de21010430c
md"
## Solution:
"

# ╔═╡ 3a809403-89e4-42c6-8dd3-f9dd6c915d14
reactions = [
    [-1, 0, 0, 0, 0, 0],
    [0, -1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, -1, -1, 1, 1, -1],
    [-1, 0, 1, -1, -1, 1],
    [0, 1, 1, -1, -1, 1],
    [1, 0, -1, 1, 1, -1],
]

# ╔═╡ af84b4f4-ae44-4e9c-95a0-0d591d425ba1
function calculate_propensity(population::Vector{Int}, rates::Vector{Float64})
    propensity = Vector{Int}(undef, 8)
    propensity[1:4] = population[1:4] .* rates[1]
    propensity[3:4] = population[3:4] .* rates[2]
    propensity[5] = population[3] * population[2] * rates[3]
    propensity[6] = population[4] * population[1] * rates[3]
    propensity[7] = population[5] * rates[4]
    propensity[8] = population[6] * rates[4]
    return propensity
end

# ╔═╡ f06290ce-ce15-41a7-8700-37653627753a
function simulate(
    init_population::Vector{Int},
    rates::Vector{Float64},
    reactions::Vector{Vector{Int}},
    N::Int,
)
    population = Array{Int}(undef, length(init_population), N)
    population[:, 1] = init_population
    for i = 2:N
        propensity = calculate_propensity(population[:, i-1], rates)
        d = Categorical(propensity ./ sum(propensity))
        population[:, i] = population[:, i-1] .+ reactions[rand(d)]
    end
    return population
end

# ╔═╡ 0d4c4807-9821-422e-8ea7-b7618efbc060
init_population_a = [0, 0, 1, 0, 0, 1]

# ╔═╡ 49d581a2-9c0d-42d3-8bed-bb27b34088e5
rates = [2.0, 10.0, 2.0, 1.0]

# ╔═╡ e366458e-f5c1-44bf-b188-750307370ada
begin
    population = simulate(init_population_a, rates, reactions, 1000)
    hline([0], line = (4, :dash, :darkgreen), lab = "")
    plot!(
        1:1000,
        population[1, :] .- population[2, :],
        line = (:steppre, :steelblue),
        lab = L"n_A-n_B",
    )
end

# ╔═╡ Cell order:
# ╠═1398e536-09ff-11ec-177f-6b6668e95c95
# ╠═21524504-7968-4a8d-b99f-1de21010430c
# ╠═3a809403-89e4-42c6-8dd3-f9dd6c915d14
# ╠═af84b4f4-ae44-4e9c-95a0-0d591d425ba1
# ╠═f06290ce-ce15-41a7-8700-37653627753a
# ╠═0d4c4807-9821-422e-8ea7-b7618efbc060
# ╠═49d581a2-9c0d-42d3-8bed-bb27b34088e5
# ╠═e366458e-f5c1-44bf-b188-750307370ada
