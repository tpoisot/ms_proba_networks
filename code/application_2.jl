using Gadfly
using DataFrames

include("proba_utils.jl")
include("matrix_utils.jl")
include("degree.jl")
include("connectance.jl")
include("katz.jl")
include("nestedness.jl")

## Read the Robertson data

function null1(A)
   return ones(A) .* connectance(A)
end

function null2(A)
   ki = degree_out(A) ./ size(A)[2]
   kj = degree_in(A) ./ size(A)[1]
   return (ki .+ kj')./2.0
end

A = readdlm("../data/robertson.txt")
t1 = null1(A)
t2 = null2(A)
