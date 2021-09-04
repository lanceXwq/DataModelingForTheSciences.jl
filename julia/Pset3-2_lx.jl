### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 1398e536-09ff-11ec-177f-6b6668e95c95
begin
    import Pkg, Random
    Pkg.activate("")
    Random.seed!(1234)
    using Distributions, Plots, LaTeXStrings
end

# ╔═╡ 25d421ef-33e3-400f-ad81-f8e09d950dc7
md"
# Problem 2: The genetic toggle switch and stochastic bistability.
The toggle switch (Gardner et al. Nat. 403, 339-342) is a common feedback loop motif in systems biology and it exhibits a fascinating behavior called  'stochastic bistability'. We will now simulate this behavior. Consider the following chemical reactions involving two proteins, ``A`` and ``B``:
```math
\begin{align*}
A & \xrightarrow{k_{d}} \phi\\
B & \xrightarrow{k_{d}} \phi\\
g_{A} & \xrightarrow{k_{p}} g_{A}+A\\
g_{B} & \xrightarrow{k_{p}} g_{B}+B\\
g_{A}+B & \xleftrightarrow[k_{r}]{k_{f}} g_{A}^{*}\\
g_{B}+A & \xleftrightarrow[k_{r}]{k_{f}} g_{B}^{*}
\end{align*}
```
where ``k_{d}`` are degradation rates and ``k_{p}`` are production rates for both proteins.
``g_{A}`` is the gene responsible for the production of ``A`` which can converted into an inactive form ``g_{A}^{*}`` by binding to B.
Vice versa for ``g_{B}``. Assume you only have one gene available in the cell so that ``g_{A}+g_{A}^{*}=1`` and  ``g_{B}+g_{B}^{*}=1``.
Also, assume throughout that ``g_{A}^{*}+g_{B}^{*} =1``, ``k_{d} < k_{p}`` and ``k_{r}<k_{f}n_{B}, k_{f}n_{A}``.\
(a) Simulate the chemical reactions starting with ``n_{A}=0`` and ``n_{B}=0`` for many time steps. Adjust your rates until you see 
stochastic switching events between periods when ``A`` exceeds ``B`` in number and ``B`` exceeds ``A`` in number. 
You should see stochastic hopping between two solutions (which we call ''fixed points'').\
(b) Would you expect to see this stochastic switching occur if you had started with a large amount of ``n_{A}`` and ``n_{B}`` initially? 
In technical language, qualitatively explain (in words) how the fixed point structure changes for the corresponding rate equations.\
(c) The condition that ``g_{A}^{*}+g_{B}^{*} =1`` is called the exclusive switch. Relax this condition and re-simulate the toggle switch.
What new fixed point appears?
"

# ╔═╡ 21524504-7968-4a8d-b99f-1de21010430c
md"
## Solution:
"

# ╔═╡ b0b548fb-f781-4e89-a5e4-bf9f337bddd1
md"
### a)
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
# ╟─25d421ef-33e3-400f-ad81-f8e09d950dc7
# ╟─21524504-7968-4a8d-b99f-1de21010430c
# ╠═1398e536-09ff-11ec-177f-6b6668e95c95
# ╟─b0b548fb-f781-4e89-a5e4-bf9f337bddd1
# ╠═3a809403-89e4-42c6-8dd3-f9dd6c915d14
# ╠═af84b4f4-ae44-4e9c-95a0-0d591d425ba1
# ╠═f06290ce-ce15-41a7-8700-37653627753a
# ╠═0d4c4807-9821-422e-8ea7-b7618efbc060
# ╠═49d581a2-9c0d-42d3-8bed-bb27b34088e5
# ╠═e366458e-f5c1-44bf-b188-750307370ada
