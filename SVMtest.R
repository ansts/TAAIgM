require(e1071)
require(plot3D)
cpl2=colorRampPalette(c("#AFAFAF0A","#FF0FFF0A","#00FFF00A"))

X=vn10
colnames(X)=paste(dgn,1:21, sep="")
prX=sapply(1:21,function(i){        
  y=cmdscale(dist(t(X[FSbootsC[[i]],])))   
  s=svm(y[-i,], as.factor(dgn[-i]), gamma = 0.07, cost=1000, probability = T)
  x=predict(s, newdata=t(y[i,]))
  r=apply(y,2,range)
  pry=sapply(seq(r[,1][1], r[,1][2], length.out = 10), function(x1){         
    sapply(seq(r[,2][1], r[,2][2], length.out = 10), function(x2){
      attr(predict(s, newdata=t(c(x1,x2)), probability = TRUE), "probabilities")[1]
    })
  })
  fnm=paste("svmG\\SVMpred",i,".pdf",sep="")
  pdf(file=fnm, width=7, height=6)
  op=par(mai=c(1,1,1,1))
  image2D(t(pry),resfac = 10, col=cm.colors(1000), xaxt="n", yaxt="n", xlab="", ylab="", clab="Probability")
  op=par(new=TRUE,mai=c(1,1,1,1.575))
  par(new=T)
  plot(y[-i,],cex=0, xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2")
  points(t(y[i,]), cex=5, xlim=r[,1], ylim=r[,2], pch=1, xlim=r[,1], ylim=r[,2])
  text(y, labels=dgn, col=dgnf, xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2")
  dev.off()
  return(x)
})

t=table(dgn,as.character(prX))
print(t)

prXp=prX=="GBM"
t=table(dgnG,prXp)
TP=t[4]; TN=t[1];FP=t[2];FN=t[3]
MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
accuracy=(TP+TN)/sum(t)
F1=2*TP/(2*TP+FN+FP)

parametersG=c(MCC=MCC, Acc=accuracy, F1=F1)

prX=sapply(1:21,function(i){        
  y=cmdscale(dist(t(X[FSbootsG[[i]],])))   
  s=svm(y[-i,], as.factor(dgn[-i]), gamma = 0.009, cost=50000, probability = T)
  x=predict(s, newdata=t(y[i,]))
  r=apply(y,2,range)
  pry=sapply(seq(r[,1][1], r[,1][2], length.out = 10), function(x1){         
    sapply(seq(r[,2][1], r[,2][2], length.out = 10), function(x2){
      attr(predict(s, newdata=t(c(x1,x2)), probability = TRUE), "probabilities")[1]
    })
  })
  fnm=paste("svmC\\SVMpred",i,".pdf",sep="")
  pdf(file=fnm, width=7, height=6)
  op=par(mai=c(1,1,1,1))
  image2D(t(pry),resfac = 10, col=cm.colors(1000), xaxt="n", yaxt="n", xlab="", ylab="", clab="Probability")
  op=par(new=TRUE,mai=c(1,1,1,1.575))
  par(new=T)
  plot(y[-i,],cex=0, xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2")
  points(t(y[i,]), cex=5, xlim=r[,1], ylim=r[,2], pch=1, xlim=r[,1], ylim=r[,2])
  text(y, labels=dgn, col=dgnf, xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2")
  dev.off()
  return(x)
})

t=table(dgn,as.character(prX))
print(t)

prXp=prX=="cntr"
t=table(dgnC,prXp)
TP=t[4]; TN=t[1];FP=t[2];FN=t[3]
MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
accuracy=(TP+TN)/sum(t)
F1=2*TP/(2*TP+FN+FP)
parametersC=c(MCC=MCC, Acc=accuracy, F1=F1)
