Loo=function(part, l=zsets3calls){
require(parallel)
x=seq_along(l)[((part-1)*11+1):(part*10+1)]
cl=makeCluster(16)
clusterExport(cl, c("rfe", "BHgamma", "bhgamfix","connfix","dunnfix","vb10","dgnC","dgnG","dgnM","dgnf","l","x"), envir=environment())
fs=parSapply(cl,x, function(i){
  m=vb10[l[[i]],-i]
  calls=lapply(list(dgnC, dgnM, dgnG), function(d){
      d=d[-i]
      x=rfe(m, d)
      rownames(x[[1]])
  })
  x=table(unlist(calls))
  names(x)[x>1]
})
stopCluster(cl)
save(fs, file=paste(c("F",part), collapse="_"))
return(fs)
}
