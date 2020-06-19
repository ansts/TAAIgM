# See descripton in scanconn()

critsel=function(D,s){
  require(Rfast)
  D=rowsums(s^D)
  return(c(which.max(D),max(D)))
}