#! /usr/bin/env julia

using ProbabilisticNetwork

# Read the Poullain data
A = readdlm("../data/poullain.txt")

# Draw them in a PNG file
# draw_matrix(A; file="../figures/poullain.png")

# 1000 Bernoulli trials
distr_links = [A |> make_bernoulli |> links for i in 1:1000]
distr_nest = [A |> make_bernoulli |> nestedness for i in 1:1000]

# Average nestedness
mapslices(mean, hcat(distr_nest...), 2)

# Motif 1
b1 = make_unipartite([1.0 0.0; 0.0 1.0])
b2 = make_unipartite([1.0 1.0; 0.0 1.0])
b3 = make_unipartite([1.0 1.0; 1.0 1.0])

uA = make_unipartite(A)
println("PROBA")
println("b1:", motif(uA, b1), " b2:", motif(uA, b2), " b3:", motif(uA, b3))

uA = make_binary(A)
println("BINARY")
println("b1:", motif(uA, b1), " b2:", motif(uA, b2), " b3:", motif(uA, b3))


#println(motif(make_unipartite(A), m1))

# Bernoulli trials
#nm1 = map((x) -> motif(make_unipartite(x), m1), nullmodel(A, n=50))
#nm1prime = map((x) -> motif(make_unipartite(x), m1'), nullmodel(A, n=50))
