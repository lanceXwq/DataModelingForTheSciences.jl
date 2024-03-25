### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ a8969360-c690-11ec-16de-15eb0f10d48c
begin
    import Pkg
    Pkg.activate(Base.current_project())
    using GLMakie
    using Random
    #using Distributions
    using LinearAlgebra
    Random.seed!(1234)
end

# ╔═╡ 9c6eee60-6506-432f-84c9-ef6e1a0a0932
function emission(data, μ, σ)
    return [exp(-(data - μ[1])^2 / (2 * σ^2)), exp(-(data - μ[2])^2 / (2 * σ^2))]
end

# ╔═╡ 49003c95-ced7-4284-81e5-af949b4a8a17
π = [0.97 0.03; 0.05 0.95]

# ╔═╡ 172b8b3a-2e9f-46ea-beb0-da1823fa7e19
π_transposed = transpose(π)

# ╔═╡ c2bacc8b-7d89-44f7-8105-dbccbfe19a10
ρ = [0.5; 0.5]

# ╔═╡ dc7c6df7-ba19-4267-b7ef-ee511b73b641
N = 300

# ╔═╡ 74d174e0-ddc4-4de8-868a-5f1ec3a62a25
μ = [1, 2]

# ╔═╡ e6126127-b641-4fc8-bb76-182b0968c235
σ = 0.2

# ╔═╡ 7282a265-7b29-40fb-b7fc-784aecf1ea6d
s = Vector{Int8}(undef, N)

# ╔═╡ a3a19d14-d9fe-4269-a21f-d13485a45e8e
begin
    for i = 2:N
        s[i] = s[i-1] - (rand() > π[s[i-1], s[i-1]])
        s[i] == 0 && (s[i] = 2)
    end
end

# ╔═╡ 4b614a2a-2883-4e8a-8db7-b5b9a2519d9c
data = μ[s] .+ σ .* randn(length(s))

# ╔═╡ 305bce1f-5c9e-4e1a-9ba8-952acac46376
begin
    stairs(s, color=:black)
    scatter!(data)
    current_figure()
end

# ╔═╡ 48c0efff-18e6-403d-b3e1-4a1e3757ec92
A = Array{Float64,2}(undef, 2, N)

# ╔═╡ 0675add7-fafe-484f-8edc-6e2b591e76a1
begin
    A[:, 1] = emission(data[1], μ, σ)
    A[:, 1] = A[:, 1] ./ sum(A[:, 1])
    for i = 2:N
        A[:, i] = emission(data[i], μ, σ) .* (π_transposed * view(A, :, i - 1))
        A[:, i] = A[:, i] ./ sum(A[:, i])
    end
end

# ╔═╡ d0f6e22e-20c1-429a-8efc-6b50dd92ec6a
viterbi_path = Vector{Int8}(undef, N)

# ╔═╡ 3c37d0a5-99bc-4ebf-bacd-193a061969cf
begin
    viterbi_path[N] = 1 + (A[2, N] > A[1, N])
    for i = N-1:-1:1
        viterbi_path[i] = 1 + (A[2, i] * π[2, viterbi_path[i+1]] > A[1, i] * π[1, viterbi_path[i+1]])
    end
end

# ╔═╡ 3df423ea-07d2-4e0c-9aa4-fdc33584b774
begin
    stairs!(viterbi_path, color=:red)
    current_figure()
end

# ╔═╡ Cell order:
# ╠═a8969360-c690-11ec-16de-15eb0f10d48c
# ╠═9c6eee60-6506-432f-84c9-ef6e1a0a0932
# ╠═49003c95-ced7-4284-81e5-af949b4a8a17
# ╠═172b8b3a-2e9f-46ea-beb0-da1823fa7e19
# ╠═c2bacc8b-7d89-44f7-8105-dbccbfe19a10
# ╠═dc7c6df7-ba19-4267-b7ef-ee511b73b641
# ╠═74d174e0-ddc4-4de8-868a-5f1ec3a62a25
# ╠═e6126127-b641-4fc8-bb76-182b0968c235
# ╠═7282a265-7b29-40fb-b7fc-784aecf1ea6d
# ╠═a3a19d14-d9fe-4269-a21f-d13485a45e8e
# ╠═4b614a2a-2883-4e8a-8db7-b5b9a2519d9c
# ╠═305bce1f-5c9e-4e1a-9ba8-952acac46376
# ╠═48c0efff-18e6-403d-b3e1-4a1e3757ec92
# ╠═0675add7-fafe-484f-8edc-6e2b591e76a1
# ╠═d0f6e22e-20c1-429a-8efc-6b50dd92ec6a
# ╠═3c37d0a5-99bc-4ebf-bacd-193a061969cf
# ╠═3df423ea-07d2-4e0c-9aa4-fdc33584b774
