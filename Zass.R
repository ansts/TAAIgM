Zass=function(G, membrshp){
  require(igraph)
  X=assortativity(G, membrshp, weights=edge.attributes(G)$weight)
  r=sapply(1:1000, function(i) {
    assortativity(G, sample(membrshp), weights=edge.attributes(G)$weight)
  })
  c(X,(X-mean(r))/sd(r)) # assortativity and z-score
}