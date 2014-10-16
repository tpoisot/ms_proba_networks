using Gadfly
using DataFrames
using PEN

searchdir(path,key) = filter(x->contains(x,key), readdir(path))

data_folder = "../data/pollination/"

data_files = map((x) -> data_folder * x, searchdir(data_folder, ".txt"))

A_list = map(make_binary, map(readdlm, data_files))
A_1= map(null1, A_list)
A_2= map(null2l, A_list)

## Get connectance for all
Co = map(connectance, A_list)

