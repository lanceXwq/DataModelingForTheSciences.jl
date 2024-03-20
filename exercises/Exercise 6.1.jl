### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# ╔═╡ 32bfca5a-bab1-11ec-3f84-9d55923fec04
begin
    import Pkg
    Pkg.activate(Base.current_project())
    using GLMakie
    using Random
    using Distributions
	using LinearAlgebra
    Random.seed!(1234)
end

# ╔═╡ 384cb4b7-1686-48c5-9957-2d5412e8d873
v = 1

# ╔═╡ 72827156-024f-45a8-9cb2-1d00ef64dc24
input_length = 21

# ╔═╡ 6c0acf6c-8e4b-409b-9a61-d7d7f3477ba9
test_length = 200

# ╔═╡ d0c9a8c3-750f-441f-ad8a-12965b30f779
target_fxn(x) = (x+2) * (x+1) * (x-0.5) * (x-1) * (x-2) 

# ╔═╡ f0e9c819-3ec4-44f1-8032-4bd197ccf626
input_x = LinRange(-2.1, 2.1, input_length)

# ╔═╡ 92dbbf6e-3f99-4417-9c38-425d2b8816b1
input_y = target_fxn.(input_x) + sqrt(v) * randn(length(input_x))

# ╔═╡ d8103500-d9c8-4095-b2c4-d062fca6561f
begin
	x = -2.2:0.01:2.2
	y = target_fxn.(x)
	lines(x, y, color = :black)
	scatter!(input_x, input_y)
	current_figure()
end

# ╔═╡ 90e217e0-c4f5-41f4-b6b6-56edcaf1219a
test_x = LinRange(-2.5, 2.5, test_length)

# ╔═╡ 5f6812f1-9ba3-43b9-9743-9971ca3ab0dc
function get_kernel(input_x::AbstractVector, test_x::AbstractVector, decay_scale::Real, corr_strength::Real)
	x = [collect(input_x); collect(test_x)]
    C = zeros(length(x), length(x))
	@simd for i = 1:length(x)
    	C[1:i, i] = x[i] .- view(x,1:i)
        C[1:i, i] = .- view(C, 1:i,i) .^ 2 ./ decay_scale ^ 2
        C[1:i, i] = exp.(view(C, 1:i,i)) .* corr_strength ^ 2        
    end
    return(Symmetric(C))
end

# ╔═╡ f01a84f5-4495-4b89-b436-a867839c203b
C = get_kernel(input_x,test_x, 1, 1)

# ╔═╡ 088c4618-7b3b-4b8a-8d0c-f745a4be7f5a
view(C,1:input_length, 1:input_length) + v * I

# ╔═╡ 3dd87070-653d-45e6-855e-1b2a9dc8a76c
inv_kernel = inv(Symmetric(view(C,1:input_length, 1:input_length)) + v * I)

# ╔═╡ 46cfa2cc-0831-4bc0-919e-608ac6c0b2bb
C̃ = (@view C[input_length+1:end, input_length+1:end]) - (@view C[input_length+1:end,1:input_length]) * inv_kernel * (@view C[1:input_length, input_length+1:end])

# ╔═╡ fa0e9d09-313d-4587-92a5-30da03f0fbd1
maximum(C̃ - Symmetric(C̃))

# ╔═╡ d64a297c-0787-4132-96a8-547fa816eb1d
issymmetric(C̃)

# ╔═╡ a28cf1c4-5c5d-474f-abcb-02b144bf1383
μ̃ = (@view C[input_length+1:end,1:input_length]) * inv_kernel * input_y

# ╔═╡ a6f86e2a-d841-4e84-b8f7-396fc24ac0a7
GP = MultivariateNormal(μ̃, Symmetric(C̃ + 1e-14I))

# ╔═╡ e91af3f7-4f88-4df2-b300-5fdc6bf6e655
samples = rand(GP, 300)

# ╔═╡ 2befb2ac-afdd-4e27-97ed-da845bfc88d0
begin
	for sample in eachcol(samples)
		lines!(test_x, sample, color = (:red, 0.05))
	end
	current_figure()
end

# ╔═╡ fb030027-3c4a-45de-9a28-da4b42a55adc
current_axis().

# ╔═╡ 4c0e4aba-317e-43b2-a957-e59b9dd71beb


# ╔═╡ Cell order:
# ╠═32bfca5a-bab1-11ec-3f84-9d55923fec04
# ╠═384cb4b7-1686-48c5-9957-2d5412e8d873
# ╠═72827156-024f-45a8-9cb2-1d00ef64dc24
# ╠═6c0acf6c-8e4b-409b-9a61-d7d7f3477ba9
# ╠═d0c9a8c3-750f-441f-ad8a-12965b30f779
# ╠═f0e9c819-3ec4-44f1-8032-4bd197ccf626
# ╠═92dbbf6e-3f99-4417-9c38-425d2b8816b1
# ╠═d8103500-d9c8-4095-b2c4-d062fca6561f
# ╠═90e217e0-c4f5-41f4-b6b6-56edcaf1219a
# ╠═5f6812f1-9ba3-43b9-9743-9971ca3ab0dc
# ╠═f01a84f5-4495-4b89-b436-a867839c203b
# ╠═088c4618-7b3b-4b8a-8d0c-f745a4be7f5a
# ╠═3dd87070-653d-45e6-855e-1b2a9dc8a76c
# ╠═46cfa2cc-0831-4bc0-919e-608ac6c0b2bb
# ╠═fa0e9d09-313d-4587-92a5-30da03f0fbd1
# ╠═d64a297c-0787-4132-96a8-547fa816eb1d
# ╠═a28cf1c4-5c5d-474f-abcb-02b144bf1383
# ╠═a6f86e2a-d841-4e84-b8f7-396fc24ac0a7
# ╠═e91af3f7-4f88-4df2-b300-5fdc6bf6e655
# ╠═2befb2ac-afdd-4e27-97ed-da845bfc88d0
# ╠═fb030027-3c4a-45de-9a28-da4b42a55adc
# ╠═4c0e4aba-317e-43b2-a957-e59b9dd71beb
