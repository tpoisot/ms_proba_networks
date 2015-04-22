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
