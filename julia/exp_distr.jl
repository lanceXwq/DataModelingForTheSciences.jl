# From Example 1.13 Simulating from an exponential distribution
#
# 		p(r)   = λexp(-λr) if r >= 0, 0 otherwise
# 		C(r)   = 1 - exp(-λ r)
# 		C-1(u) = -1/λ log(1 - u)
# Generate u from a uniform distibution 0 and 1
# Compute  r = -1/λ log ( 1 - u)
#
using Random, PyPlot, LinearAlgebra

function myFit(x::Array, y::Array)
	M = [sum(y) sum(@. x*y); sum(@. x*y) sum(@. x^2 * y)]
	b = [sum(@. y*log(y)); sum(@. x*y*log(y))]
	param = inv(M)*b
	return exp(param[1]), param[2]
end



function myMean(x::Array)
	N = length(x)
	mean = 0
	for el in x
		mean += el
	end
	return mean/N
end

function generateData(λ::Float64, size::Int)
	u = rand(size)
	return @. (-1/λ) .* log(1 .- u)
end


function expDistr(x, λ::Float64, type::String)
	if type == "pdf"
		return @. λ .* exp(-λ .* x)
	elseif type == "cdf"
		return @. 1 .- exp(-λ .* x)
	else
		println("Not an option!")
	end
end


r = generateData(1.0, 100)
r_sorted = sort(r) # for cdf
p = (1:length(r)) ./ (length(r) - 1)
p = convert(Array, p)
pdf_exact = expDistr(LinRange(0, 5, 100), 1.0, "pdf")
cdf_exact = expDistr(LinRange(0, 5, 100), 1.0, "cdf")

fig, axs = subplots(1, 2)
axs[1].hist(r, bins=5, density=true, label = "data")
axs[1].plot(LinRange(0, 5, 100), pdf_exact, lw = 1, label = "exact")
axs[1].set_ylabel("PDF")
axs[2].plot(r_sorted, p, label = "data")
axs[2].plot(LinRange(0, 5, 100), cdf_exact, lw = 1, label = "exact")
axs[2].set_ylabel("CDF")
legend(loc = "best")
show()

x = convert(Array, LinRange(minimum(r), maximum(r), length(r)))


A, B = myFit(x, r)


using Printf

@printf "Estimating λ...\n"
@printf "| Method | Estimate | Time |\n"
@printf "|--------------------------|\n"
@printf "| Mean   | %0.2f    |      |\n" myMean(r)
@printf "| Fit    | %0.2f    | 	    |\n" A
