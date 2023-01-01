Zmod=function(G, membrshp){
  require(igraph)
  X=modularity(G, membrshp, weights=edge.attributes(G)$weight)
  r=sapply(1:1000, function(i) modularity(G, sample(membrshp), weights=edge.attributes(G)$weight))
  c(X,(X-mean(r))/sd(r)) # modularity and z-score
}