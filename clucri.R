
# A funciton which measures the cluster separation using a composite clustering criterion 
# consisting of the best value selected among Dunn's or BHgamma criteria as well as negative Connectivity.
# The criteria are transformed to compensate for their dependence on the dimensionality 
# and to restrict them in the range between -1 and 1. The conversion functions dunnfix, connfix and bhgamfix 
# are vectorized for faster calculations.

clucri=function(m0,c){
  require(clValid)
  require(cluster)
  require(Rfast)

  if (is.logical(c)) cc=as.integer(c)+1 else cc=as.integer(c)
 
  m=m0
  N0=nrow(m0)
  nm=rownames(m0)
  nc=unique(c)
  vconf=Vectorize(connfix)
  vbhgf=Vectorize(bhgamfix)
  vdunf=Vectorize(dunnfix)
  N=nrow(m)
  d1=Dist(t(m))
  x=c(BHgamma(d1,cc), dunn(as.dist(d1),cc),connectivity(as.dist(d1),cc))
  D=sum(vbhgf(x[1]),vdunf(x[2],N),-vconf(x[3]))

  return(D)
}
