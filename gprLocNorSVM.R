#
# Reads fin gpr file with p(ath) and writes 
# the fout with B(ackground) columns equal 
# to the calculated local baseline. Flags: wr(ite), sh(o)w 
# and R(ed/Green)sw(ap) - to make the green channel 
# under red column names for pepStat. The baseline is calculated 
# using LOESS  
#
gprLocNorSVM<-function(fin, p=NULL, fout=NULL, wr=T, shw=F, Rsw=F) {
  require(parallel)
  require(reshape2)
  require(stringr)
  require(e1071)
  finp=paste(p,fin,sep = "")
  fcon=file(finp)
  fhead=readLines(con = fcon, n=34)
  f=read.delim(finp, skip=34, header=T, check.names = F)
  close(fcon)
  dfim=data.frame(cbind(f$`Row`, f$`Column`, f$`F532 Median`))
  colnames(dfim)=c("R","C","V")
  N=nrow(dfim)
  nr=max(dfim[,1])
  nc= max(dfim[,2])
  if (shw==T){
    img0=acast(dfim, R~C, value.var = "V")
    image(img0, main=c(fin," Original"), zlim=c(0,65000), col=topo.colors(128))
    lr=readline()
  }
  
  rs=7
  cs=7
  p=3
  co=1 
  for (i in 1:N){ 
    if (dfim[i,1] %in% rs:(nr-rs) & dfim[i,2] %in% cs:(nc-cs)) {                                                                  # For each spot take the spots in a patch (rs.2+1)Rx(cs.2+1)C around it
      z=dfim[abs(dfim[,1]-dfim[i,1])<=rs & abs(dfim[,2]-dfim[i,2])<=cs,]  # and keep only the bottom up to the 1/pth part .
      z1=z[rank(z[,3])<(nrow(z)/p),]                                      # Put the in btm. 
      if (co==1) {btm=z1} else {btm=rbind(btm,z1)}                        #
      btm=unique(btm)
      co=co+1
    }  
  }

  lmd=svm(btm[,1:2], btm[,3], cost = 3000, gamma=1, epsilon = 0.0001)                                     

  prd=predict(lmd,dfim[,1:2])                                       
  lv=min(dfim[,3]-prd)-1                                               
  pred=prd+lv
  lnew=as.double(dfim[,3]-pred)
  
  #ldfim=log2(lnew)
  #lmdl=svm(dfim[,1:2], ldfim, cost = 3000, gamma=1, epsilon = 0.0001)
  #my=median(ldfim)
  #yy=predict(lmdl,dfim[,1:2])
  #dfimnew=2^(ldfim-yy+my)#+pred

  if (shw==T){
    for (i in 0:3) {
      x=c((1+i*(N%/%4)),((1+i)*(N%/%4)))
      plot(dfim[x[1]:x[2],3], cex=.3, ylim = c(0, 65000))
      par(new=T)
      plot(pred[x[1]:x[2]], cex=.1,col=2, ylim = c(0, 65000))
      lines(pred[x[1]:x[2]], cex=.3,col=2, ylim = c(0, 65000))
      par(new=F)
      lr=readline()
      plot(log2(lnew[x[1]:x[2]]), cex=.4)
      lr=readline()
#      plot(log2(dfim[x[1]:x[2],3]), cex=.4)
#      lr=readline()
#      plot(dfimnew[x[1]:x[2]], cex=.4)
#      lr=readline()
      }
    hist(lnew,breaks=300)
    lr=readline()
    bkg=matrix(pred, dim(img0), byrow = T)
    image(bkg, main=c(fin," Background"), zlim=c(0,65000), col=topo.colors(128))
    lr=readline()
    #img1=img0-bkg
    img1=matrix(lnew, dim(img0), byrow = T)
    image(img1, main=c(fin," Original"), zlim=c(0,65000), col=topo.colors(128))
    
    #img1[img1<0]=0
    #image(img1, main=c(fin," Filtered"), zlim=c(0,65000), col=topo.colors(128))
    print(max(img1))
    lr=readline()
  }
  if (wr) {
    f$`B532`=pred
    f$`B532 Median`=pred
    f$`B532 Mean`=pred

    #f$`F532`=dfimnew
    #f$`F532 Median`=dfimnew
    #f$`F532 Mean`=dfimnew

    if (Rsw==T) {
      fcn=colnames(f)
      fcn=str_replace_all(fcn, "635","630")
      fcn=str_replace_all(fcn, "532", "635")
      colnames(f)=fcn
    }
    newp=paste("LoessBkgr/",sep = "")                       # Write the new .grp files in 
    dir.create(newp)                                        # subfolder LoessBkgr/.
    fout=paste(newp,fin,sep = "")
    fconw=file(fout, 'w')
    writeLines(fhead,con=fconw)
    write.table(f,file=fconw, row.names=F, sep = '\t')
    close(fconw)
  }
}