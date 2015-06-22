#! /usr/bin/env julia

using DataFrames
using ProbabilisticNetwork

searchdir(path,key) = filter(x->contains(x,key), readdir(path))

nest = (x) -> nestedness(x)[1]
mod = (x) -> x |> make_unipartite |> modularity |> Q

data_folder = "../data/canaria/"

proba_files = map((x) -> data_folder * x, searchdir(data_folder, "-proba.txt"))
binar_files = map((x) -> data_folder * x, searchdir(data_folder, "-binary.txt"))

W = map(readdlm, proba_files)
A = map(readdlm, binar_files)

conn_fix = map(connectance, A)
conn_pro = map(connectance, W)

nest_fix = map(nest, A)
nest_pro = map(nest, W)

mod_fix = map(mod, A)
mod_pro = map(mod, W)

df = DataFrame(i=1:length(W),
    bC = conn_fix, pC = conn_pro, eC = conn_fix .- conn_pro,
    bN = nest_fix, pN = nest_pro, eN = nest_fix .- nest_pro,
    bM = mod_fix, pM = mod_pro, eM = mod_fix .- mod_pro)
writetable("../figures/app3.dat", df, separator = '\t', header = true)
