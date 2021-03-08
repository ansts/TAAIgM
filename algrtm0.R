# "Pipeline"
#
# Initial raw data processing
#
# filepr preprocesses .gpr files usin gprLocNorSVM to locally 
#   normalize them by SVR, the fitted bachkgrounf is recorded 
#   as B635 columns.
# cn holds the controls' sequences
# Use pepStat to get pSet and pnSet (normlized for compositional bias in binding).
#   No spots are flagged -
# 
#
# v0[N spots x cases] is data from pnSet, 
# v1[N peptides x cases] - aggregated data, control peptides are removed  
#   as well as 12 peptides with CV exceeding 20% for at least one measurement
# vn[sames v2] - v1 after global normalization with cyclic loess 
#   from limma, method "affy".
# vb - vn after batch compesnsation
#  
# pepts - peptide list from pSet (all duplicates, same as in the chip layout)
# pepf - peptides as factor used to preserve order as in the chip as well as 
#    as in the proteins scanned
# pep - final set of sequences with order preserved 
# dgn - diagnosis

rm(list=ls())
ls()
require(stringi)
require(matrixStats)
require(pepStat)
require(limma)
require(parallel)
require(sva)
require(mixtools)
require(reshape2)
require(gplots)
require(ggsci)
require(umap)
require(clusterCrit)
require(rgl)
require(abind)
require(plot3D)
require(corrplot)
require(factoextra)
require(stringdist)
require(qualV)
require(igraph)
require(factoextra)
require(pbapply)
require(Rfast)
require(boot)
require(multimode)

cpl=colorRampPalette(c("#0000A000","#0000FF00","#00FF0000","#FFFF0000","#FFFFFF00","#FF000000"))
cpl1=colorRampPalette(c("#0000AF00","#FF000000"))

filepr("IgM/")      #
mapFile="mapall.csv"
dirToParse="SVMBkgr"
cn=c("YPYDVPDYAG", "DYKDDDDKAS")
pSet=makePeptideSet(path=dirToParse, mapping.file = mapFile)   # also applies normexp backgorund correcion based on the background calculated by gprLocNorSVM
pnSet=normalizeArray(pSet, centered = FALSE)                   # removes dependence on the chemical composition of the aa independednt of the sequence
coln=paste(1:21,pnSet@phenoData@data$diag, sep="_")
dgn=pnSet@phenoData@data$diag
dgnf=as.double(as.factor(dgn))
dgnG=dgn=="GBM"
pepts=pnSet@featureRange@elementMetadata@listData$peptide
pepf=factor(pepts, levels = unique(pepts))
v0=pnSet@assayData$exprs
v_agg_mean=aggregate(v0, by=list(pepts=pepf), mean)
v_agg_sd=aggregate(v0, by=list(pepts=pepf), sd)
v_cv=as.matrix(v_agg_sd[,-1]/v_agg_mean[,-1])
badrows=which(v_cv>0.2, arr.ind = T)[,1]
v1=as.matrix(v_agg_mean[,-1])
rownames(v1)=v_agg_mean[,1]
colnames(v1)=coln
cnrn=which(rownames(v1) %in% cn)

v1=v1[-union(cnrn, badrows),]    # remove control peptides and peptides with CV>0.2 in at least 1 pat.
pep=rownames(v1)

vn=normalizeCyclicLoess(v1, iterations=1, method="affy")
batch=c(rep(1,10),rep(2,11))
mm=model.matrix(~as.factor(dgn))
vb10=ComBat(vn, batch, mod=mm)

pdf(file="vb10dotplots.pdf", width=10, height=10)
rng=range(vb10)
plot(as.data.frame(vb10),cex=0.1, pch=16, xlim=rng, ylim=rng)
dev.off()

# Biological information about the peptides - frames encompassing two different 
# sequences are removed (they are the result of concatenating all protein sequences
# and extracting the 15-mer frames shifted by 2 residues in one pass). The file
# pep_info.csv contains the original data provided by the manufacturer.

pepinfo=read.csv(file="pep_info0.csv", header = F, row.names = NULL)  
pepinfo=unique(pepinfo)
pepinfotab=table(pepinfo[,1])
pepi=names(pepinfotab)
pepinf=pepinfo[,1]
names(pepinf)=pepinfo[,2]
pepinf=pepinf[rownames(vb10)]

pepiused=pepinf[!is.na(pepinf)]
vb10ag=vb10[names(pepiused),]
vb10agL=melt(vb10ag)
colnames(vb10agL)=c("Seq_","D_","Val")
x=stri_extract_all(vb10agL$D, regex="(?<=_)\\w+" )
vb10agL$D_=unlist(x)
x=pepiused[vb10agL$Seq_]
vb10agL=cbind(vb10agL,Ag_=x)
glm_vb10=glm(Val~D_*Ag_, family = "gaussian", data=as.data.frame(vb10agL))
sumglm=summary(glm_vb10)

sink("glm_summary.txt")
x=sumglm$coefficients[sumglm$coefficients[,4]<0.05,]
print(x[order(abs(x[,3]), decreasing = T),])
sink()

write.csv(vb10, file="vb10.csv")

plot(rowMeans(vb10),rowSds(vb10), cex=0.3)

callslm=rownames(sumglm$coefficients)[sumglm$coefficients[,4]<0.05]
callslm=callslm[grep("\\bD",callslm)]
callslm=stri_extract_all(callslm,regex="(?<=:Ag_).+")
pcallslm=pepinfo[pepinfo$V1 %in% callslm,2]
vb10s=vb10[rownames(vb10) %in% pcallslm,]  

###############################################################################
# Search for relevant profiles using z scores
###############################################################################

N=nrow(vb10)
proct=proc.time()
cl <- makeCluster(6)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cl, ex)
clusterEvalQ(cl, library(matrixStats))
clusterExport(cl, list("vb10","dgnf"), envir=environment())

zsets3=pbsapply(1:21, function(i){             # leave one out loop
  X=vb10[,-i]
  dg=as.double(dgnf[-i])
  dgu=unique(dg)

  y=sapply(1:20,function(a){ 
    dgi=dgu[dgu!=dg[a]]
    pm=sapply(dgi, function(dia){
      i=dg==dia
      m=scale(t(cbind(X[,a],X[,i])))[1,]
    })
    rowMaxs(pm, value = T)
  })
  rowSums(abs(y>2))},  cl=cl )                                      
stopCluster(cl)                                  
print(proc.time()-proct)                         

rownames(zsets3)=rownames(vb10)
zsets3calls=apply(zsets3,2,function(col) rownames(zsets3)[col>1])

save(vb10,dgnG, dgnf, zsets3calls, file="vars")

#zsets3=array(zsets3,dim = c(N,20,21))

#zsets3calls=sapply(1:21,function(i){
#  m=zsets3[,,i]
#  x=rowSums(m>.45)==3
#  rownames(vb10)[x]
#})

confusiomM_direct=SVMforValidation(zsets3calls)
for (i in 1:21) {plotMDS(vb10[zsets3calls[[i]],], col=dgnf, main=coln[i])}

proct=proc.time()
cl <- makeCluster(6)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cl, ex)
clusterExport(cl, list("vb10","dgnf","zsets3calls"), envir=environment())
zsets3callsfs=pblapply(1:21, function(i){   # Recursive feature selection step
  l=zsets3calls[[i]]
  x=rfe(vb10[l,-i], dgnf[-i])
  return(rownames(x[[1]]))
}, cl=cl)
stopCluster(cl)
print(proc.time()-proct) 

load("F_1")
zsets3callsfs=fs
load("F_2")
zsets3callsfs=c(zsets3callsfs,fs)

confusiomM_fs=SVMforValidation(zsets3callsfs)


# The profiles are ordered by the number of "bootstrap" sets the features  
# are common for (using the leave one out scheme data as a bootstrap of a sort) 
# starting from none (all features are concatenated) to max 21 
# - the features found in all (or most) bootstrap sets. 

FSzsetst=table(unlist(zsets3callsfs))
FSzscom=sapply(0:20, function(i){
  x=names(FSzsetst)[FSzsets3t>i]
  x[order(FSzsetst[x])]
})

for (i in 1:21) plotMDS(vb10[FSzscom[[i]],], col=dgnf, main=i)

FSzscomcri=sapply(1:20, function(i){  # Selection of the bootstrap set 
  clucri(vb10[FSzscom[[i]],], dgnf)
})

pdf(file="hmzSets.pdf")
heatmap.2(vb10[FSzscom[[5]],], hclustfun = hclwrd, scale="row", na.rm=F,key.title = NA, colsep=c(2,5,10,15),rowsep=c(10,15,21,31,42,52,60,82),symkey=FALSE, cexCol=1,cexRow=0.4, trace="none",col=cpl1(100), margins = c(5,12), lwid=c(0.5,1), lhei = c(0.2,1))
dev.off()  #  

hclFSzsc=hclust(dist(vb10[FSzscom[[5]],]), method="ward.D2")
ct15=cutree(hclFSzsc,h=15)
ct10=cutree(hclFSzsc,h=10)
ct6=cutree(hclFSzsc,h=6)
cts=cbind(ct15,ct10,ct6)
cts=cts[hclFSzsc$order,]
pepFS5=rownames(cts)
pats=colnames(vb10)

clstspep=c()
clstspep[[1]]=rownames(cts)[cts[,1]==4]
clstspep[[2]]=rownames(cts)[cts[,1]==1]
clstspep[[3]]=rownames(cts)[cts[,2]==5]
clstspep[[4]]=rownames(cts)[cts[,2]==3]
clstspep[[5]]=rownames(cts)[cts[,2]==6]
clstspep[[6]]=rownames(cts)[cts[,3]==4]
clstspep[[7]]=rownames(cts)[cts[,3]==12]
clstspep[[8]]=rownames(cts)[cts[,3]==21]
clstspep[[9]]=rownames(cts)[cts[,3]==5]

clstprots=lapply(clstspep,function(i) pepiused[i])
GBMup=sort(unlist(clstprots[c(5,8,9)]))
GBMdown=sort(unlist(clstprots[c(2,4,6)]))

unqprots=table(pepiused)
x=sapply(unqprots,function(x) x=0)
gbmupt=table(GBMup)
x[names(gbmupt)]=gbmupt
gbmupt=x
x=sapply(unqprots,function(x) x=0)
gbmdnt=table(GBMdown)
x[names(gbmdnt)]=gbmdnt
gbmdnt=x
ttltbl=cbind(unqprots-(gbmdnt+gbmupt),gbmupt,gbmdnt)
xsqttltb=chisq.test(ttltbl, simulate.p.value = T)
sort(xsqttltb$residuals[,2])
sort(xsqttltb$residuals[,3])

z=ttltbl[rowSums(ttltbl[,2:3])>0,]
