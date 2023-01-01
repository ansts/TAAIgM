Zass_pro=function(G=Gvs, mx=vs10, prof){
  require(igraph)
  require(parallel)
  require(pbapply)
  
  mx=mx-min(mx)+1
  vG=names(V(G))
  mx=mx[vG,]
  crit=sapply(vG, function(p){
    pp=names(ego(G,1,p)[[1]])
    x=mx[pp,]
    print(which(vG==p))
    if (!is.null(nrow(x))) {
      n=nrow(x)
      Y=colMeans(x)
    }
    else {
      n=1
      Y=x
    }
    xA=as.factor(prof)
    diff(aggregate(Y, by=list(xA),"mean")[,2])
  })
  
  X=assortativity(G, crit)
  cl=makeCluster(15)
  clusterExport(cl,c("G","vG","mx","prof"))
  clusterEvalQ(cl,require(igraph))
  r=pbsapply(1:150, function(i) {
    prof=sample(prof)
    crit=t(sapply(vG, function(p){
      pp=names(ego(G,1,p)[[1]])
      x=mx[pp,]
      print(which(vG==p))
      if (!is.null(nrow(x))) {
        n=nrow(x)
        Y=colMeans(x)
      }
      else {
        n=1
        Y=x
      }
      xA=as.factor(prof)
      diff(aggregate(Y, by=list(xA),"mean")[,2])
    }))
    assortativity(G, crit)
  }, cl=cl)
  stopCluster(cl)
  c(X,(X-mean(r))/sd(r)) # assortativity and z-score
}