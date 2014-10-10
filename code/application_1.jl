using Gadfly
using DataFrames

include("proba_utils.jl")
include("matrix_utils.jl")
include("degree.jl")
include("connectance.jl")
include("katz.jl")
include("nestedness.jl")

## Read the Poullain data

A = make_unipartite(readdlm("../data/poullain.txt"))

distr_links = float([A |> make_binary |> float |> links for i in 1:1000])
distr_nest = float([A |> make_binary |> float |> nestedness for i in 1:1000])
