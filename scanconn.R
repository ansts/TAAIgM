
# A funciton which iteratively scans the set of features and finds at each cycle the feature  
# the removal of which imrpoves most the separation between the classes of interest c, then removes it 
# leaving one less features for the next cycle. The separation is measured using a composite clustering criterion in critsel()
# consisting of the best value selected among Dunn's or BHgamma criteria as well as negative Connectivity.
# The criteria are transformed to compensate for their dependence on the dimensionality in the range of dimensions 3-600 
# and to restrict them in the range between -1 and 1. The conversion functions dunnfix, connfix and bhgamfix 
# are vectorized for faster calculations. The dimesionality standardized values obtained are converted using exponential function 
# with a base s randomly selected in the range of 1.5 - 5 and common  for all calculations in one call to this function.
# This transformation ensures the non-negativity of the criteria and adds a small random factor selecting in every call
# a slightly different profile. In this way the huge profile space is searched more efficiently in a kind of bootstrap 
# fashion. The 3 criteria thus transformed are summed up to produce the composite criterion. 
# The final profile is determined by consensus of several (n=5) profiles in profileSearch().

scanconn=function(m0,c, s=NULL){
  require(clValid)
  require(parallel)
  require(cluster)
  require(Rfast)
  require(reshape2)
  require(matrixStats)
  require(beepr)
  if (is.logical(c)) cc=as.integer(c)+1 else cc=as.integer(c)
  #pt=proc.time()
  D0=c()
  Ds0=c()
  bad=list()
  bads=list()
  if (is.null(s)) s=runif(1,1.5,5) 
  print(s)
  m=m0
  nm=rownames(m0)
  nc=unique(c)
  vconf=Vectorize(connfix)
  vbhgf=Vectorize(bhgamfix)
  vdunf=Vectorize(dunnfix)
  repeat {
      D=c()
      N=nrow(m)
      print(N)
      if (N<3) break
      D=t(sapply(1:N, function(i) {
            d1=Dist(t(m[-i,]))
            x=c(BHgamma(d1,cc), dunn(d1,cc),connectivity(d1,cc))
            return(x)
        }))
      NN=rep(N-1,N)
      D=cbind(vbhgf(D[,1],NN),vdunf(D[,2],NN),-vconf(D[,3],NN))
      D=critsel(D,s)
      md=D[1]                    
      dmin=D[2]
      bad=c(bad,rownames(m)[md])
      m=m[-md,]
      D0=rbind(D0,c(N,dmin)) #,quantile(D[,2],0.000125)
  }
  
  bad=unlist(bad)
  yr=range(D0[,2])
  if (nrow(D0)>1) plot(nrow(D0):1,D0[,2])
  if (length(D0)>0) crit=which.max(D0[,2]) else return(list(NA,NA,NA))
  #print(crit)
  if (crit<2) crit=2
  m=m0[!(nm %in% bad[1:(crit-1)]),]
  
  return(list(m,bad,D0))
}