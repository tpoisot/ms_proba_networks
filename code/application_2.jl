#! /usr/bin/env julia

using DataFrames
using ProbabilisticNetwork

searchdir(path,key) = filter(x->contains(x,key), readdir(path))

prepare = (x) -> x |> make_binary
nest = (x) -> nestedness(x)[1]

data_folder = "../data/pollination/"

data_files = map((x) -> data_folder * x, searchdir(data_folder, ".txt"))

A_list = map(prepare, map(readdlm, data_files))
A_1 = map(null1, A_list)
A_2 = map(null2, A_list)
A_3i = map(null3in, A_list)
A_3o = map(null3out, A_list)

## Get nestedness for all
NO = map(nest, A_list)

dNO1 = NO .- map(nest, A_1)
dNO2 = NO .- map(nest, A_2)
dNO3i = NO .- map(nest, A_3i)
dNO3o = NO .- map(nest, A_3o)

df = DataFrame(i=1:length(data_files), N=NO, d1 = dNO1, d2 = dNO2, d3i = dNO3i, d3o = dNO3o)
writetable("../figures/app2.dat", df, separator = '\t', header = true)
