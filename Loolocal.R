Loolocal0=function(l,m, oo, d1=dgnC, d2=dgnM, d3=dgnG){
  require(parallel)
  require(pbapply)
 
  L=list(d1,d2,d3)
  cl=makeCluster(3)
  clusterExport(cl,list("BHgamma","bhgamfix","connfix","dunnfix","rfe","L", "m","l","oo"), envir=environment())
  calls=pbsapply(L, function(d) {
    m0=m[l,-oo]
    x=rfe(m0, d[-oo])
    rownames(x[[1]])
  },cl=cl)
  x=table(unlist(calls))
  fs=names(x)[x>1]
  stopCluster(cl)
  return(fs)
}
