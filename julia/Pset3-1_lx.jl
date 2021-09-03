### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 97e973bc-5903-434f-b1b2-238c7ddcdf71
begin
    import Pkg, Random
    Pkg.activate("")
	Random.seed!(1234)
    using Distributions, Plots
end

# ╔═╡ dc618f4e-b069-11eb-0755-db6734f11bef
md"
# Problem 1: Stochastic binary decisions.
Cells fate decisions are often based on small initial fluctuations that are amplified and reinforced (through feedback) over time. 
These events are called stochastic binary decisions (Artyomov et al., PNAS, 104, 18958, 2007). 
We will now simulate such a process. Consider the following set of chemical rxns:
```math
\begin{align*}
A_1&\xrightarrow{k_1}E+A_1\\
A_2&\xrightarrow{k_2}S+A_2\\
E+A_1&\xrightarrow{k_3}E+A_1^*\\
A_1^*&\xrightarrow{k_4}E+A_1^*\\
S+A_1&\xrightarrow{k_5}S+A_1^{INACTIV}\\
S&\xrightarrow{k_d}\phi\\
E&\xrightarrow{k_d}\phi
\end{align*}
```
where ``A_1`` is an agonist and ``A_2`` is its antagonist. 
``E`` is an enzyme that converts ``A_1`` into its protected form ``A_1^*`` and ``A_1^∗ `` , in turn, stimulates the production of E (positive feedback). 
We are interested in the steady state amount of ``A_1^*``. 
If ``S`` is present it can permanently de-activate ``A_1`` .\
(a) Start with 10 agonists and 10 antagonists with ``k_1 = k_2 = k_d = k_4 = 1, k_3 = 100, k_5 = 100``.
Simulate the process to completion using Gillespie’s algorithm and histogram the final amount of ``A_1^*`` . 
Explain, in words, the result you obtain.\
(b) Repeat the simulation and histograming of part (a) starting with 1000 agonists and 1000 antagonists.
Explain, in words, how your histogram differs from that of part (a).\
(c) If you had solved the corresponding rates equations explain in words what you would expect the steady state population of ``A_1^*`` to look like.
"

# ╔═╡ 3a8ebc91-5596-481b-8438-817462899448
md"
## Solution:
"

# ╔═╡ b003b560-e394-47aa-956e-72dbb1ca5927
reactions = [
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [-1, 0, 0, 0, 1, 0],
    [0, 0, 1, 0, 0, 0],
    [-1, 0, 0, 0, 0, 1],
    [0, 0, -1, 0, 0, 0],
    [0, 0, 0, -1, 0, 0],
]

# ╔═╡ 2220fdb3-72b8-46f9-b685-c3cceb3ac427
rates = [1., 1., 100., 1., 100., 1.]

# ╔═╡ 16c0951f-eb96-49e7-91cb-da99e04c0fcc
function calculate_propensity(population::Vector{Int}, rates::Vector{Float64})
    propensity = zeros(7)
    propensity[1:2] = population[1:2] .* rates[1:2]
    propensity[3] = population[3] * population[1] * rates[3]
    propensity[4] = population[5] * rates[4]
    propensity[5] = population[4] * population[1] * rates[5]
    propensity[6] = population[3] * rates[6]
    propensity[7] = population[4] * rates[6]
    return propensity
end

# ╔═╡ b4d3faa9-4dc7-470b-82b9-e7dc97b45bb5
function simulate(init_population::Vector{Int}, rates::Vector{Float64}, rxns::Vector{Vector{Int}})
    population = [init_population]
    while population[end][1] > 0
        propensity = calculate_propensity(population[end], rates)
        d = Categorical(propensity ./ sum(propensity))
        push!(population, population[end])
        population[end] += rxns[rand(d)]
    end
    return population[end][5]
end

# ╔═╡ 5cf00f2c-f4a1-4aa3-8afe-93697cfdc4a2
md"
### a)
"

# ╔═╡ fd393590-58dd-4a25-867e-74f487b7a4d3
init_population_a = [10, 10, 0, 0, 0, 0]

# ╔═╡ 45f34728-fb4c-4ccd-a815-76ff1a682642
let
    final_population = zeros(Int, 2000)
    for idx = 1:length(final_population)
        final_population[idx] = simulate(init_population_a, rates, reactions)
    end
    histogram(final_population, norm = true, label = "", bin = 10, xticks = 0:11)
end

# ╔═╡ af96a06e-1a95-4cd9-aef1-78a01858f376
md"
### b)
"

# ╔═╡ cc87baa1-93e7-4a54-af3b-e48e7d1bbce2
init_population_b = [1000, 1000, 0, 0, 0, 0]

# ╔═╡ 94bcd1dc-7bc9-401c-b912-591f1f071b38
let
    final_population = zeros(Integer, 2000)
    for idx = 1:length(final_population)
        final_population[idx] = simulate(init_population_b, rates, reactions)
    end
    histogram(final_population, norm = true, label = "")
end

# ╔═╡ Cell order:
# ╟─dc618f4e-b069-11eb-0755-db6734f11bef
# ╟─3a8ebc91-5596-481b-8438-817462899448
# ╠═97e973bc-5903-434f-b1b2-238c7ddcdf71
# ╟─b003b560-e394-47aa-956e-72dbb1ca5927
# ╟─2220fdb3-72b8-46f9-b685-c3cceb3ac427
# ╠═16c0951f-eb96-49e7-91cb-da99e04c0fcc
# ╠═b4d3faa9-4dc7-470b-82b9-e7dc97b45bb5
# ╟─5cf00f2c-f4a1-4aa3-8afe-93697cfdc4a2
# ╟─fd393590-58dd-4a25-867e-74f487b7a4d3
# ╠═45f34728-fb4c-4ccd-a815-76ff1a682642
# ╟─af96a06e-1a95-4cd9-aef1-78a01858f376
# ╠═cc87baa1-93e7-4a54-af3b-e48e7d1bbce2
# ╠═94bcd1dc-7bc9-401c-b912-591f1f071b38
