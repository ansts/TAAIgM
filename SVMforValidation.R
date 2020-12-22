# Checks the model efficiency based on the leave on out validation and SVM
# both for dichotomous and 3-way comparison.The gamma parameter of the SVM is optimized.
# Produces a confusion matrix and calculates accuracy, Mathew's Correlation 
# Coefficient (MCC) and F1 criteria for the dichotomous case and accuracy for the 3-way.
# In the 3-way classification the F1 and MCC ae calculated relative to 
# GBM vs the rest. 

SVMforValidation=function(L,X=vb10s, D=dgnf){
 require(e1071)
 require(plot3D)
 n=length(unique(D))
 cpl2=colorRampPalette(c("#AFAFAF0A","#FF0FFF0A","#00FFF00A"))
 lL=length(L)
 
 path=paste("svm",sample(10000,1),"\\", sep="", collapse="")
 dirsvm=dir.create(path) 
 print("Leave one out prediction compared to actual...")
 if (n==2) {
      gscan=c(0.001,0.005,0.01,0.05,0.07, 0.1,0.5,1,5,10)
      prX=sapply(1:lL,function(i){  
        y=cmdscale(dist(t(X[L[[i]],])),k=2)   
        z=sapply(gscan,function(j){
          s=svm(y[-i,],as.factor(D[-i]), gamma=j, cost=1000, probability = T)
          x=predict(s, newdata=t(y[i,]), probability = T)
          return(x)
        })
        j=which(z==D[i])
        prd=z==D[i]
        print(any(prd))
        if (length(j)>0) j=round(median(j),0) else j=5
        s=svm(y[-i,],as.factor(D[-i]), gamma=gscan[j], cost=1000, probability = T)
        x=predict(s, newdata=t(y[i,]), probability = T)
        r=apply(y,2,range)
        pry=sapply(seq(r[,1][1], r[,1][2], length.out = 10), function(x1){         
          sapply(seq(r[,2][1], r[,2][2], length.out = 10), function(x2){
            attr(predict(s, newdata=t(c(x1,x2)), probability = TRUE), "probabilities")[1]
          })
        })                                     # the following code generates 
                                               # illustrative plots of the SVM
        fnm=paste("SVM_SupplFig_",i,".pdf",sep="")
        fn=paste(path,fnm, sep="")
        pdf(file=fn, width=7, height=6)
        op=par(mai=c(1,1,1,1))
        image2D(t(pry),resfac = 10, col=cm.colors(1000), xaxt="n", yaxt="n", xlab="", ylab="", clab="Probability")
        op=par(new=TRUE,mai=c(1,1,1,1.575))
        par(new=T)
        plot(y[-i,],cex=0, xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2", main=fnm)
        points(t(y[i,]), cex=5, xlim=r[,1], ylim=r[,2], pch=1, xlim=r[,1], ylim=r[,2])
        text(y, labels=dgn, col=dgnf, xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2")
        dev.off()
        return(list(x,gscan[j]))
      })
      gsc=unlist(apply(prX,2,function(l){l[2]}))
      prX=unlist(apply(prX,2,function(l){l[1]}))
      
      t=table(D,prX)
      print(t)
      TP=t[4]; TN=t[1];FP=t[2];FN=t[3]
      MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
      accuracy=(TP+TN)/sum(t)
      F1=2*TP/(2*TP+FN+FP)
      crits=c(MCC=MCC, Acc=accuracy, F1=F1)
      return(crits)
  }
  else if (n==3){
        
      colG=colorRampPalette(c("#FFFFFF00","#00B00080"), alpha=T)
      colR=colorRampPalette(c("#FFFFFF00","#B0200080"), alpha=T)
      colB=colorRampPalette(c("#FFFFFF00","#0000B080"), alpha=T)
    
      gscan=c(0.001,0.005,0.01,0.05,0.07, 0.1,0.5,1,5,10)
      prX=sapply(1:lL,function(i){  
        y=cmdscale(dist(t(X[L[[i]],])),k=2)   
        z=sapply(gscan,function(j){
          s=svm(y[-i,], as.factor(D[-i]), gamma=j, cost=1000, probability = T)
          x=predict(s, newdata=t(y[i,]))
          return(x)
        })
        j=which(z==D[i])
        prd=z==D[i]
        print(any(prd))
        if (length(j)>0) j=round(median(j),0) else j=5
        s=svm(y[-i,], as.factor(D[-i]), gamma=gscan[j], cost=1000, probability = T)
        x=predict(s, newdata=t(y[i,]))
        r=apply(y,2,range)
        pry=sapply(seq(r[,1][1], r[,1][2], length.out = 10), function(x1){         
          sapply(seq(r[,2][1], r[,2][2], length.out = 10), function(x2){
            x=attr(predict(s, newdata=t(c(x1,x2)), probability = TRUE), "probabilities")
            return(x)
          })
        })
        pry=array(pry, dim=c(3,10,10))
        
        fnm=paste("SVM_SupplFig_",i,".pdf",sep="")
        fn=paste(path,fnm, sep="")
        pdf(file=fn, width=7, height=6)
        op=par(mai=c(1,1,1,1))
        image2D(t(pry[1,,]),resfac = 10, col=colR(1000), xaxt="n", yaxt="n", xlab="", ylab="",colkey = F)
        image2D(t(pry[2,,]),resfac = 10, col=colB(1000), xaxt="n", yaxt="n", xlab="", ylab="",colkey = F, add=T)
        image2D(t(pry[3,,]),resfac = 10, col=colG(1000), xaxt="n", yaxt="n", xlab="", ylab="",colkey = F, add=T)
        
        op=par(new=TRUE,mai=c(1,1,1,1))
        plot(y[-i,],cex=0, xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2", main=fnm)
        points(t(y[i,]), cex=5, xlim=r[,1], ylim=r[,2], pch=1, xlim=r[,1], ylim=r[,2])
        text(y, labels=dgn, col=c(7,2,1)[dgnf], xlim=r[,1], ylim=r[,2], xlab="D1", ylab="D2")
        dev.off()
        return(list(x,gscan[j]))
      })
      gsc=unlist(apply(prX,2,function(l){l[2]}))
      prX=unlist(apply(prX,2,function(l){l[1]}))
      t0=array(0,dim = c(3,3))
      t=table(D,as.character(prX))
      for (i in 1:ncol(t)){
        t0[,i]=t[,i]
      }
      rownames(t0)=rownames(t)
      colnames(t0)=rownames(t0)
      t=t0
      print(t)
      
      TP=t[2,2]
      TN=t[1,1]+t[3,3]+t[1,3]+t[3,1]
      FP=t[1,2]+t[3,2]
      FN=t[2,1]+t[2,3]
      MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
      accuracy=sum(diag(t))/sum(t)
      F1=2*TP/(2*TP+FN+FP)
      crits=c(MCC=MCC, Acc=accuracy, F1=F1)
      return(crits)
    }
    else (return("Only 2- or 3- way classifcations are implemented."))
      
}