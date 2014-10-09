#=

Makes a bipartite network unipartite

=#
function make_unipartite(A::Array{Float64,2})
   S = sum(size(A))
   B = zeros(Float64,(S,S))
   B[1:size(A)[1],size(A)[1]+1:S] = A
   return B
end
