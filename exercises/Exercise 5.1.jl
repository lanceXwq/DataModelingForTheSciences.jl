### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ 7a51cb04-9e41-11ec-2c7c-39e816fd07c3
begin
	import Pkg
	Pkg.activate(Base.current_project())
	using GLMakie, Random
	Random.seed!(1234)
end

# ╔═╡ 7a6605b9-07bf-445c-b6ee-f916304c215f
r = rand(1000, 3) .* 2 .- 1

# ╔═╡ 7c96fadf-9555-44ee-b17f-d50e574493c0
if_real = @. r[:,2] ^ 2 - 4 * r[:,1] * r[:,3] > 0

# ╔═╡ 3ce3a293-c9c3-4a75-a04c-778f2b0c4ed5
set_theme!(backgroundcolor = :white)

# ╔═╡ f3003a97-3d4d-4b36-ae29-515f961adac2
scatter(r[:,1],r[:,2],r[:,3], color = if_real, colormap = :jet, markersize = 20)

# ╔═╡ aebc94a5-a5b5-4c77-87e0-93dfaa1c4c4a
count(if_real) / 1000

# ╔═╡ Cell order:
# ╟─7a51cb04-9e41-11ec-2c7c-39e816fd07c3
# ╠═7a6605b9-07bf-445c-b6ee-f916304c215f
# ╠═7c96fadf-9555-44ee-b17f-d50e574493c0
# ╠═3ce3a293-c9c3-4a75-a04c-778f2b0c4ed5
# ╠═f3003a97-3d4d-4b36-ae29-515f961adac2
# ╠═aebc94a5-a5b5-4c77-87e0-93dfaa1c4c4a
