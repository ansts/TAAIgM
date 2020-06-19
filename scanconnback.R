# A function for forward selection of features starting typically from a set selected by scanconn() followed by scanning 
# the set of remaining features in the original data matrix. At each cycle the feature which imrpoves most the separation 
# between the classes of interest c, is added to the selected feature set and is removed from the feature set to further scan. 
# The separation is measured as described in the scanconn(). To shorten the calculations a criterion based 
# on the slope of the last third of the curve of the criterion dependence on feature set length is used to interrupt 
# the calculations after a maximum has been reached. The feature set corresponding to the maximal criterion is returned.

scanconnback=function(m0,c,pl,s=NULL){
  require(clValid)
  require(parallel)
  require(cluster)
  require(Rfast)
  require(reshape2)
  require(matrixStats)
  require(beepr)
  if (is.logical(c)) cc=as.integer(c)+1 else cc=as.integer(c)
  pt=proc.time()
  D0=c()
  Ds0=c()
  good=list()
  if (is.null(s)) s=runif(1,1.5,5) 
  nm=rownames(m0)
  m=m0[pl,]
  m1=m0[!(nm %in% pl),]
  n=nrow(m)
  nc=unique(c)
  fail=F
  vconf=Vectorize(connfix)
  vdunf=Vectorize(dunnfix)
  vbhgf=Vectorize(bhgamfix)
  #x0c=connectivity(Dist(t(m)),cc)
  #x0d=dunn(dist(t(m)),cc)
  nD=0
  repeat {
      D=c()
      N=nrow(m1)
      print(N)
      n=nrow(m)
      if (length(D0)>6) {
        nD=nrow(D0)
        crit=line(D0[(nD%/%3):nD,])
        crit=crit$coefficients[2]
        }
      else crit=0
      if (N<3|(crit<(-0.001)&(nD>45))) break  
      D=t(sapply(1:N, function(i) {
            m2=rbind(m,m1[i,])  
            d1=Dist(t(m2))
            x=c(BHgamma(d1,cc), dunn(d1,cc),connectivity(d1,cc))
            return(x)
        }))
      nn=rep(n+1,N)
      D=cbind(vbhgf(D[,1],nn),vdunf(D[,2],nn),-vconf(D[,3],nn))
      D=critsel(D,s)
      md=D[1]                    
      dm=D[2]
      d0=dm
      good=c(good,rownames(m1)[md])
      m=rbind(m,m1[md,])
      rownames(m)[nrow(m)]=rownames(m1)[md]
      m1=m1[-md,]
      D0=rbind(D0,c(n,dm)) 
  }
  
  if (max(D0[,2]) <= d0) return (list(m[pl,],NULL, NULL)) 
  
  good=unlist(good)
  yr=range(D0[,2])
  m=m[1:(length(pl)+which.max(D0[,2])),]
  
  return(list(m,good,D0))
}