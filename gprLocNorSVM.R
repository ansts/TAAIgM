#
# Reads fin gpr file with p(ath) and writes (if wr is TRUE)
# the fout with B(ackground) columns equal 
# to the calculated local baseline. Flags: wr(ite), sh(o)w 
# and R(ed/Green)sw(ap) - to put the green channel data
# under the red column names for pepStat to use (pepStat takes the red column). 
# The baseline is calculated using support vector regression.
# Show flag shw allows for plotting intermediate data.
#

gprLocNorSVM<-function(fin, p=NULL, fout=NULL, wr=T, shw=F, Rsw=F) {
  require(parallel)
  require(reshape2)
  require(stringr)
  require(e1071)
  
  finp=paste(p,fin,sep = "")
  #fcon=file(finp)
  f2l=readLines(finp, n=2)
  nlns=as.double(str_extract(f2l[2], "[1-9]+"))
  fhead=readLines(finp, n=nlns+2)
  f=read.delim(finp, skip=nlns+2, header=T, check.names = F)
  #close(fcon)
  dfim=data.frame(cbind(f$`Row`, f$`Column`, f$`F532 Median`))
  colnames(dfim)=c("R","C","V")
  N=nrow(dfim)
  nr=max(dfim[,1])
  nc=max(dfim[,2])
  if (shw==T){
    img0=acast(dfim, R~C, value.var = "V")
    image(img0, main=c(fin," Original"), zlim=c(0,65000), col=topo.colors(128))

  }
  
  rs=7
  cs=7
  p=3
  co=1 
  for (i in 1:N){ 
    if (dfim[i,1] %in% rs:(nr-rs) & dfim[i,2] %in% cs:(nc-cs)) {          # For each spot take the spots in a patch of 
      z=dfim[abs(dfim[,1]-dfim[i,1])<=rs & abs(dfim[,2]-dfim[i,2])<=cs,]  # (rs x 2 + 1) x (cs x 2 + 1) spots around it
      z1=z[rank(z[,3])<(nrow(z)/p),]                                      # and keep the those bellow 1st 1/p-quantile.
      if (co==1) {btm=z1} else {btm=rbind(btm,z1)}                        # Put them in btm. 
      btm=unique(btm)
      co=co+1
    }  
  }

  lmd=svm(btm[,1:2], btm[,3], cost = 3000, gamma=1, epsilon = 0.0001)    # cost and gamma have been otpimized                                  

  prd=predict(lmd,dfim[,1:2])                                            
  lv=min(dfim[,3]-prd)-1                                                 
  pred=prd+lv                                                            # extracted baseline (background)
  lnew=as.double(dfim[,3]-pred)

  if (shw==T){
    for (i in 0:3) {
      x=c((1+i*(N%/%4)),((1+i)*(N%/%4)))
      plot(dfim[x[1]:x[2],3], cex=.3, ylim = c(0, 65000))
      par(new=T)
      plot(pred[x[1]:x[2]], cex=.1,col=2, ylim = c(0, 65000),  ylab="")
      lines(pred[x[1]:x[2]], cex=.3,col=2, ylim = c(0, 65000),  ylab="")
      par(new=F)
      plot(log2(lnew[x[1]:x[2]]), cex=.4,  ylab="")
    }
    hist(lnew,breaks=300)
    bkg=matrix(pred, dim(img0), byrow = T)
    image(bkg, main=c(fin," Background"), zlim=c(0,65000), col=topo.colors(128))
    img1=matrix(lnew, dim(img0), byrow = T)
    image(img1, main=c(fin," Original"), zlim=c(0,65000), col=topo.colors(128))
    print(paste("Max intensity: ",round(max(img1),0),sep="",collapse=""))
  }
   if (wr) {
    f$`B532`=pred
    f$`B532 Median`=pred
    f$`B532 Mean`=pred

    if (Rsw==T) {
      fcn=colnames(f)
      fcn=str_replace_all(fcn, "635","630")
      fcn=str_replace_all(fcn, "532", "635")
      colnames(f)=fcn
    }
    newp=paste("SVMBkgr/",sep = "")                       # Write the new .grp files in 
    dir.create(newp)                                        # subfolder SVMBkgr/.
    fout=paste(newp,fin,sep = "")
    fconw=file(fout, 'w')
    writeLines(fhead,con=fconw)
    write.table(f,file=fconw, row.names=F, sep = '\t')
    close(fconw)
    return()
  }
}