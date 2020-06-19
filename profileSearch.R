# See description in scanconn().

profileSearch=function(M,sff,sfb,code,nm,cy=5){
  require(rgl)
  require(parallel)
  cl <- makeCluster(5)
  ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
  clusterEvalQ(cl, {library(clValid)})
  clusterExport(cl, ex)
  clusterExport(cl, list("M","sff","sfb", "code"), envir=environment())
  
  xi=parLapply(cl,1:cy, function(t) {
    prct=proc.time()
    x=scanconn(M[sff,], code)
    y1=rownames(x[[1]])
    xb=scanconnback(M[sfb,],code,y1)
    y2=rownames(xb[[1]])
    print(proc.time()-prct)
    return(list(y1,y2))
  })
  
  stopCluster(cl)
  x=list()
  x[[1]]=lapply(xi,function(l){l[[1]]})
  x[[2]]=lapply(xi,function(l){l[[2]]})
  x=lapply(x, function(l){table(unlist(l))})
  
  x5=lapply(x, function(l){names(l[l==5])})
  xl=scanconn(M[x5[[2]],], code)
  if (length(xl[[2]])<4) return(x5[[2]]) else return(rownames(xl[[1]]))
}