# A wrapper for the limma linear model based extraction of calls

limfit1_<-function(fdsn,vals, n=1000, p=0.05, ip=FALSE, cmf=1){
  require(limma)
  require(sva)
  require(parallel)
  m=length(unique(fdsn))
  ncl=ncol(vals)
  if (ncl!=length(fdsn)) {return("Error: ncol of design and matrix differ!")}
  clin=character(0)
  design = model.matrix(~0+fdsn)
  cn=paste("f",0:(m-1),sep ="")
  colnames(design)=cn
  cntr=sapply(1:(m-1),function(i){paste("f",i,"-f0", sep="")})
  cmx=makeContrasts(contrasts=cntr, levels=cn)
  fit=lmFit(vals,design)
  fit2=contrasts.fit(fit,cmx)
  fit2=eBayes(fit2, trend = F)
  #fit1=eBayes(fit)
  tptb=topTable(fit2, coef=c(1,2), number=n) #,adjust.method="holm"
  nr=nrow(tptb[tptb$adj.P.Val<p,])
  clin=rownames(tptb[1:nr,])
  res=list(clin,tptb)
  return(res)
}  

