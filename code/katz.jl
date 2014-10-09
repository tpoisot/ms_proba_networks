function katz_centrality(A::Array{Float64,2}; a::Float64=0.1, k::Int64=5)
   @assert a <= 1.0
   @assert a >= 0.0
   @assert k >= 1
   centr = sum(hcat(map((k) -> vec(sum((a^k).*(A^k),1)), [1:5])...),1)
   return centr ./ sum(centr)
end
