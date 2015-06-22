library(betalink)
library(rmangal)
library(plyr)

api = mangalapi()

dataset = getDataset(api, 1)

A = alply(dataset$network, 1, function(x) toIgraph(api, x), .progress=create_progress_bar(name="text"))

M = metaweb(A)

plants = which(degree(M, mode="in")>0)
pollinators = which(degree(M, mode="out")>0)

pairs = data.frame(expand.grid(plant=names(plants), pollinator=names(pollinators)))
pairs$presence = rep(0, nrow(pairs))
pairs$interaction = rep(0, nrow(pairs))

for(i in c(1:nrow(pairs))) {
  pair = pairs[i,]
  for(a in A) {
    sp = V(a)$name
    if((pair$plant %in% sp) & (pair$pollinator %in% sp))
    {
      pair$presence = pair$presence + 1
      al = data.frame(get.edgelist(a))
      if(nrow(subset(al, (X1 == as.vector(pair$pollinator)) & (X2 == as.vector(pair$plant)))))
      {
        pair$interaction = pair$interaction + 1
      }
    }
  }
  pairs[i,] = pair
}

pairs$probability = pairs$interaction / pairs$presence
pairs$probability[is.nan(pairs$probability)] = 0

m = xtabs(probability~plant+pollinator, pairs)

write.table(m, file="../data/canaria/metaweb.txt", row.names=F, col.names=F)

for(i in c(1:length(A))) {
  a = t(get.adjacency(A[[i]], sparse=F))
  a[a>0] = 1
  plants = which(degree(A[[i]], mode="in")>0)
  pollinators = which(degree(A[[i]], mode="out")>0)
  w = m[which(rownames(m) %in% names(plants)), which(colnames(m) %in% names(pollinators))]
  a = a[which(rownames(a) %in% names(plants)), which(colnames(a) %in% names(pollinators))]
  write.table(w, file=paste("../data/canaria/", i, "-proba.txt", sep=''), row.names=F, col.names=F)
  write.table(a, file=paste("../data/canaria/", i, "-binary.txt", sep=''), row.names=F, col.names=F)
}
