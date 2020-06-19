# Analysis of teh nomalized intensities in vn1

require(stringi)
require(matrixStats)
require(pepStat)
require(limma)
require(prallel)
require(sva)
require(mixtools)
require(reshape2)
require(d3heatmap)
require(gplots)
require(ggsci)
require(umap)
require(clusterCrit)
require(rgl)
require(Rfast)
require(plot3D)
require(corrplot)
require(factoextra)
require(stringdist)
require(qualV)
require(igraph)

load("vn1")
pepinfo0=read.csv(file="pep_info0.csv", header = F, row.names = NULL) # pep_info0.csv contains the biological information about the peptides
pepinfo=unique(pepinfo0)
pepinfotab=table(pepinfo[,1])
pepi=names(pepinfotab)

cpl=colorRampPalette(c("#0000A000","#0000FF00","#00FF0000","#FFFF0000","#FFFFFF00","#FF000000"))
cpl1=colorRampPalette(c("#0000AF00","#FF000000"))


dgn=c("GBM","GBM","GBM","GBM","GBM","meta","meta","cntr","cntr","cntr","GBM","GBM","GBM","GBM","meta","meta","cntr","cntr","cntr","meta","GBM")
dgnf=as.factor(dgn)
dgnf=as.double(dgnf)
dgnG=dgn=="GBM"
dgnC=dgn=="cntr"
dgnM=dgn=="meta"

# determining the negative set by clustering using umap to visualize cluster with minimal intensity in all patients 

uconf=umap.defaults
uconf$n_neighbors=5
umvn10=umap(vn1, config = uconf)
plot(umvn10$layout,col=cpl1(4524)[rank(rowMeans(vn10))])
negi=umvn10$layout[,1]>1 & umvn10$layout[,2]<(-6)
plot(umvn10$layout[negi,],col=cpl1(4524)[rank(rowMeans(vn10))], xlim=range(umvn10$layout[,1]),ylim=range(umvn10$layout[,2]))
vn0=mean(vn1[negi,])
vn10=vn1-vn0

pepinf=pepinfo[,1]
names(pepinf)=pepinfo[,2]
pepinf=pepinf[rownames(vn10)]

x0=c(mean(vn10[negi,]),sd(vn10[negi,])) # means and SDs of the negative set of peptides (by peptide)
vbz=sapply(1:21, function(j){           # z-scores of the peptides relative to the negative peptides 
  z0=(vn10[,j]-x0[1])/(x0[2]/sqrt(length(vn10[negi,])))
  return(z0)
})

vbz0p=t(apply(vbz,1,function(i){    # transform to p
  pnorm(i,lower.tail=F)
}))

vbz0p=p.adjust(vbz0p)               # false discovery adjustment
vbz0p=matrix(vbz0p,nrow = 4524)
z0=apply(vbz0p,1,function(p){       # take only peptides (z0 is a boolean flag) with more than 10 patients                                               significant at level p
  sum(p<0.001)>10
})

# z scores of each peptide reactivtiy in each patient relative to the mean intensities in the two alternative diagnostic groups of patients for the same peptide 

ndg=table(dgnf)
x=t(vn10)
xm=aggregate(x,by=list(dgnf), FUN="mean")
xm=t(xm)[-1,]
xsd=aggregate(x,by=list(dgnf), FUN="sd")
xsd=t(xsd)[-1,]
x=t(x)

y=sapply(1:21,function(c){
  z=sapply(1:3, function(j){
    (x[,c]-xm[,j])/(xsd[,j]*sqrt(ndg[j]))
  })
  z=z[,-dgnf[c]]
  return(rowMax(abs(z)))     # take only the higher absolute value of the two z-scores 
})

z1=rowMeans(y)  # for each peptide take the mean of the maximal z-scores for the different patients

bf=rownames(vn10[z0&(z1>0.85),]) 
bb=rownames(vn10[z0&(z1>0.67),])

# Chi square test for the preference for the selected peptide set bf in the peptides from the different antigens
 
bfinfo=sapply(bf,function(p){pepinfo[pepinfo[,2]==p,1]})
bfinfo=bfinfo[lengths(bfinfo)==1]
nbf=rownames(vn10)[!(rownames(vn10) %in% bf)]
nbfinfo=sapply(nbf,function(p){pepinfo[pepinfo[,2]==p,1]})
bfinftab=t(sapply(pepi, function(p){
  Agpos=sum(bfinfo==p)
  Agneg=sum(nbfinfo==p)
  c(Agpos,Agneg)
}))
colnames(bfinftab)=c("Agpos","Agneg")
bfinftab=bfinftab[rowSums(bfinftab)>0,]
chsqbf=chisq.test(bfinftab, simulate.p.value = T)
pepinfotab=table(bfinfo)
stres=chsqbb$stdres[order(rowMax(abs(chsqbb$stdres)),decreasing = T),]
write.csv(stres, file="stres.csv")

# Wrapper feature selection based on clustering criteria

FSbootsG=lapply(1:21, function(i){                                    
  x=profileSearch(vn10[,-i],bf, bb, dgnG[-i])                          
  return(x)                             
})                                       

# summarize
FSbootG=table(unlist(FSbootsG))
FSbootcomG=sapply(0:20, function(i){x=names(FSbootG)[FSbootG>i]; x[order(FSbootG[x])]}) # profiles ordered by the number of bootstrap sets the features 
# are common for starting from none (all features are concatenated) to 
# 21 - the features found in all bootstrap sets (from a leave one out scheme relative to the patients)

for (i in 1:21) plotMDS(vn1[FSbootcomG[[i]],], col=dgnf+1, main=i) # inspection

# Best separation feature set - usually the set of features common for about half of the bootstrap feature sets

ij=11
theSetG=vn10[FSbootcomG[[ij]],]
theSetG=t(scale(t(theSetG),center = T))
colnames(theSetG)=paste(dgn,1:21,sep = "")
d3heatmap(theSetG, hclustfun = hclwrd)
pdf(file="heatmapG.pdf", width=6, height=12)  # Figure 2
  heatmap.2(theSetG, hclustfun = hclwrd, cexRow = 0.5, cexCol=0.6, trace="none", col=cpl1(128), key.title = NA, colsep=10, rowsep=57, margins = c(5,12), lwid=c(0.5,1), lhei = c(0.25,1) )
dev.off()
hctheSG=hclust(dist(theSetG), method = "ward.D2")
setordG=FSbootcomG[[ij]][hctheSG$order]
setinfG=sapply(setordG, function(p){pepinfo[pepinfo[,2]==p,1]})

featGtab=t(sapply(pepinfo[,2],function(p){
  pos=p %in% unlist(setordG[58:138])
  neg=p %in% unlist(setordG[1:57])
  return(c(pos,neg))
}))
featGtab=aggregate(featGtab, by=list(pepinfo[,1]), FUN=sum)
colnames(featGtab)=c("Ag", "Pos","Neg")
featGtab=cbind(featGtab,unlist(pepinfotab[featGtab$Ag]))
featGtab=featGtab[,-4]
featGtab[,4]=featGtab[,4]-(featGtab[,2]+featGtab[,3])
colnames(featGtab)[4]="Rest"
chsqG=chisq.test(featGtab[,2:4], simulate.p.value = T)
AgsGImpPos=cbind(chsqG$stdres[order(abs(chsqG$residuals[,1]),decreasing = T),1],featGtab[order(abs(chsqG$residuals[,1]), decreasing = T),1])

AgsGImpNeg=cbind(chsqG$stdres[order(abs(chsqG$residuals[,2]),decreasing = T),2],featGtab[order(abs(chsqG$residuals[,2]), decreasing = T),1])


# The same analysis for the other dichotomy - tumor bearing vs conrol

FSbootsC=lapply(1:21, function(i){                                    
  x=profileSearch(vn10[,-i],bf, bb, dgnC[-i])                          
  return(x)                             
})  
                                   
# summarize
FSbootC=table(unlist(FSbootsC))
FSbootcomC=sapply(0:20, function(i){x=names(FSbootC)[FSbootC>i]; x[order(FSbootC[x])]}) 
for (i in 1:21) plotMDS(vn10[FSbootcomC[[i]],], col=dgnf+1, main=i)

ij=12
theSetC=vn10[FSbootcomC[[ij]],]
theSetC=t(scale(t(theSetC),center = T))
theSetC=scanconn(theSetC, dgnC)[[1]]

colnames(theSetC)=paste(dgn,1:21,sep = "")
hctheSC=hclust(dist(theSetC), method = "ward.D2")
hcCi=cutree(hctheSC, k=4)
d3heatmap(theSetC, hclustfun = hclwrd, cexCol = .5)
pdf(file="heatmapClast.pdf", width=4, height=15)    #Figure 4
  heatmap.2(theSetC, hclustfun = hclwrd, cexRow = .3, cexCol=0.6, trace="none", col=cpl1(128), key.title = NA, colsep=c(5,15), rowsep=c(111,190,260), margins = c(5,6), lhei = c(.2,1), keysize = 2)
dev.off()
setordC=FSbootcomC[[ij]][hctheSC$order]
setinfC=sapply(setordC, function(p){pepinfo[pepinfo[,2]==p,1]})

featCtab=t(sapply(pepinfo[,2],function(p){
    cl1=p %in% unlist(setordC[1:80])
    cl2=p %in% unlist(setordC[81:150])
    cl3=p %in% unlist(setordC[151:229])
    cl4=p %in% unlist(setordC[230:340])
    return(c(cl4,cl3,cl2,cl1))
}))
featCtab=aggregate(featCtab, by=list(pepinfo[,1]), FUN=sum)
colnames(featCtab)=c("Ag", "Cl4","Cl3", "Cl2","Cl1")
featCtab=cbind(featCtab,unlist(pepinfotab[featCtab$Ag]))
featCtab=featCtab[,-6]
featCtab[,6]=featCtab[,6]-rowsums(as.matrix(featCtab[,2:5]))
colnames(featCtab)[6]="Rest"
featCtab2=data.frame(Cl4_3=(as.double(featCtab[,2])+as.double(featCtab[,3])),Cl2_1=(as.double(featCtab[,4])+as.double(featCtab[,5])),Rest=as.double(featCtab[,6]))
rownames(featCtab2)=featCtab[,1]
fi=apply(featCtab2,1,function(i){
  sum(as.double(i[1:3])==0)<2
})
featCtab2=featCtab2[fi,]

chsqC=chisq.test(featCtab2[,1:3], simulate.p.value = T)
AgsCImpPos=data.frame(stres=chsqC$stdres[order(abs(chsqC$residuals[,1]),decreasing = T),1],Ag=featCtab2[order(abs(chsqC$residuals[,1]), decreasing = T),1])

AgsCImpNeg=data.frame(stres=chsqC$stdres[order(abs(chsqC$residuals[,2]),decreasing = T),2],Ag=featCtab2[order(abs(chsqC$residuals[,2]), decreasing = T),1])
AgsCImpPos=AgsCImpPos[abs(AgsCImpPos[,1])>1.5,]
AgsCImpNeg=AgsCImpNeg[abs(AgsCImpNeg[,1])>1.5,]

#############

pdf(file="sepC.pdf", width = 10,height=10)  # Figure 3
plotMDS(theSetC,col=dgnf+1)
dev.off()

pdf(file="sepG.pdf", width = 10,height=10)  # Figure 1
plotMDS(theSetG,col=dgnf+1)
dev.off()

###########
###########
# Picking the most significant reactivities from the two feature sets selected so far

  f=dgnf                          
  f=factor(f)
  tSCc_=limfit1_(f,theSetC, p=0.051)
  theminsetC=vn10[tSCc_[[1]],]
  colnames(theminsetC)=colnames(theSetC)
  theminsetC=t(scale(t(theminsetC)))
  xhm=heatmap.2(theminsetC, hclustfun = hclwrd, cexRow = 0.5, cexCol=0.6, trace="none", col=cpl1(128), key.title = NA, colsep=c(6,16), rowsep=9, margins = c(5,6), lhei = c(.2,1), keysize = 2)
  tSCord=tSCc_[[1]][xhm$rowInd[19:1]]
  tSCinfo=pepinf[tSCord]
  write.csv(unique(tSCinfo), file="tSCinfo.csv")
  
  
  f=dgnf                          
  f=factor(f)
  tSGc_=limfit1_(f,theSetG, p=0.1)
  theminsetG=vn10[tSGc_[[1]],]
  colnames(theminsetG)=colnames(theSetG)
  theminsetG=t(scale(t(theminsetG)))
  xhm_g=heatmap.2(theminsetG, hclustfun = hclwrd, cexRow = 0.5, cexCol=0.6, trace="none", col=cpl1(128), key.title = NA, colsep=c(10,16), rowsep=c(13,20), margins = c(5,6), lhei = c(.2,1), keysize = 2)
  tSGord=tSGc_[[1]][xhm_g$rowInd[34:1]]
  tSGinfo=pepinf[tSGord]
  write.csv(unique(tSGinfo), file="tSGinfo.csv")
  
# Taking the union - tSGC, of the most significant features. The linkers are artificially formed sequences in the transitions from one antigen sequences to the next which are not interpretable biologically.
  
  x0=rownames(vn10)
  linkers=setdiff(x0,pepinfo[,2])
  tSGC=union(tSGord,tSCord)
  tSGC=setdiff(tSGC,linkers)
  theminsetGC=vn10[tSGC,]
  colnames(theminsetGC)=colnames(theSetG)
  theminsetGC=t(scale(t(theminsetGC)))
  plotMDS(theminsetGC, col=dgnf)
  png(file="heatmapGCminsets.png", width=6, height=8, units = "in", res=1200) # Figure 6
  xhm_gc=heatmap.2(theminsetGC, hclustfun = hclwrd, cexRow = 0.5, cexCol=0.6, trace="none", col=cpl1(128), key.title = NA, colsep=c(10,16), rowsep=c(10,20,30,41), margins = c(5,6), lhei = c(.2,1), keysize = 2)
  dev.off()
  tSGCord=tSGC[xhm_gc$rowInd[46:1]]
  tSGCinfo=pepinfo[pepinfo[,2] %in% tSGCord,]
  tSGCinfo=tSGCinfo[order(tSGCord),]
  write.csv(unique(tSGCinfo), file="tSGCinfo.csv")
  tSGCord_org=factor((tSGCord %in% tSGord)*1+(tSGCord %in% tSCord)*2, labels = c("G","C","GC"))
  
  tSGCinf=t(sapply(unique(pepinfo[,1]), function(i){
    x0=length(pepinfo[pepinfo[,1]==i,1])
    x1=length(tSGCinfo[tSGCinfo==i])
    print(c(i,(tSGCinfo[tSGCinfo==i])))
    return(c(x1,x0-x1))
  }))
  
  tSGCchsq=chisq.test(tSGCinf, simulate.p.value = T)
  tSGCchsqsr=tSGCchsq$stdres[order(abs(tSGCchsq$stdres[,1]),decreasing = T),1]
  pnorm(tSGCchsqsr, lower.tail = F)
  
  