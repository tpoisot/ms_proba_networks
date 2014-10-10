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
A = make_unipartite(A)
t1 = make_unipartite(t1)
t2 = make_unipartite(t2)


#=distr_links = float([A |> make_binary |> float |> links for i in 1:1000])=#
#=distr_nest = float([A |> make_binary |> float |> nestedness for i in 1:1000])=#
