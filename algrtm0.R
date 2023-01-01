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
require(reshape2)
require(gplots)
require(ggsci)
require(clusterCrit)
require(rgl)
require(plot3D)
require(stringdist)
require(qualV)
require(igraph)
require(factoextra)
require(pbapply)
require(boot)
require(multimode)
require(multcomp)
require(Biostrings)
require(heatmap3)
require(pROC)


# Collecting the raw data -------------------------------------------------

cpl=colorRampPalette(c("#00000080","#0000B080","#00B0B080","#00FF0A80","#FFFF0080","#FF000080"), alpha=T) #c("#00000080","#0000B080","#00B0B080","#00FF0A80","#FFFF0080","#FF000080")
cpl1=colorRampPalette(c("#10202000","#00806000","#E0C02040","#FFFF0000"))

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
dgnM=dgn=="meta"
dgnC=dgn=="cntr"
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
write.csv(vb10, file="vb10.csv")

pdf(file="vb10dotplots.pdf", width=10, height=10)
rng=range(vb10)
plot(as.data.frame(vb10),cex=0.1, pch=16, xlim=rng, ylim=rng)
dev.off()


# Biological information about the peptides - -----------------------------
# frames encompassing two different 
# sequences are removed (they are the result of concatenating all protein sequences
# and extracting the 15-mer frames shifted by 2 residues in one pass). The file
# pep_info.csv contains the original data provided by the manufacturer.

pepinfo=read.csv(file="pep_info.csv", row.names = NULL)  
pepinfo=pepinfo[,-c(1:2)]
pepinfotab=table(pepinfo[,2])
pepi=names(pepinfotab)
garbg=c("control", "0 0", "x\\d+")
garbg=unlist(sapply(garbg, function(p) grep(p,pepi)))
pepigrb=pepi[garbg]
pepi=pepi[-garbg]


# Correct inconsistencies in protein names --------------------------------

x=pepi
x[x==pepi[[23]]]=pepi[[5]]
x[x==pepi[[19]]]=pepi[[25]]
x[x %in% pepi[c(17,18,20,21,22)]]=pepi[[24]]
x[x %in% pepi[c(13,54)]]=pepi[[14]]
x[x==pepi[[37]]]=pepi[[39]]
x[x==pepi[[53]]]=pepi[[12]]

pepinfo=pepinfo[!(pepinfo$protein %in% pepigrb),]

for (i in seq_along(pepi)){
  pepinfo[pepinfo$protein==pepi[[i]],2]=x[[i]]
}

pepinfotab=table(pepinfo[,2])
pepi=names(pepinfotab)

pepinf=pepinfo[,2]
names(pepinf)=pepinfo[,1]
pepinf=pepinf[rownames(vb10)]
pepiused=pepinf[!is.na(pepinf)]

pepinfotab=table(pepiused)
write.csv(pepinfotab, file = "pepinfotab.csv")
pepinfotab_hu=read.csv(file="pepinfotab_hu.csv")
pepinfotab_hu=pepinfotab_hu[,-1]


# Reshaping the data to long form for the linear model tests --------------

vb10ag=vb10[names(pepiused),]
vb10agL=melt(vb10ag)
colnames(vb10agL)=c("Seq_","D_","Val")
x=stri_extract_all(vb10agL$D, regex="(?<=_)\\w+" )
vb10agL$D_=unlist(x)
x=pepiused[vb10agL$Seq_]
vb10agL=cbind(vb10agL,Ag_=x)
vb10agL$D_=as.factor(vb10agL$D_)
vb10agL$Ag_=as.factor(vb10agL$Ag_)
x=aggregate(vb10agL$Val, by=list(vb10agL$Ag_), "mean")
x[which.max((scale(x[,2]))),1]
Dc=rep("Tu", nrow(vb10agL))
Dc[vb10agL$D_=="cntr"]="None"
vb10agL=cbind(vb10agL, Tumor=Dc)
vb10agL=as.data.frame(vb10agL)

# LM ----------------------------------------------------------------------

glm_vb10=glm(Val~D_*Ag_, family = "gaussian", data=vb10agL)
sumglm=summary(glm_vb10)
sglm=sumglm$coefficients[sumglm$coefficients[,4]<0.05,]
j=grep("^Ag_",rownames(sglm))
sglm[j,]=sglm[j,][order(sglm[j,1],decreasing = T),]
j=(j[length(j)]+1):nrow(sglm)
sglm[j,]=sglm[j,][order(sglm[j,1],decreasing = T),]
write.csv(sglm, file="MeanReactivitiesComp.txt")

posth=summary(glht(glm_vb10, lsm(pairwise ~ D_|Ag_), by=NULL))
resposth=cbind(posth$test$coefficients,posth$test$pvalues)
colnames(resposth)=c("Coefficients","P.Value")
resposth[resposth[,2]<0.1,]

glm_vb10_Tu=glm(Val~Tumor*Ag_, family = "gaussian", data=vb10agL)
sumglm=summary(glm_vb10_Tu)
sglm=sumglm$coefficients[sumglm$coefficients[,4]<0.05,]

glm_vb10_epi_meta=glm(Val~Seq_*Tumor, family = "gaussian", data=vb10agL, subset = D_!="GBM")
sumglm_epi_meta=summary(glm_vb10_epi_meta)
pp=p.adjust(sumglm_epi_meta$coefficients[,4]) 
sumglm_epi_meta=sumglm_epi_meta$coefficients[pp<0.05,]
sumglm_epi_meta=sumglm_epi_meta[order(sumglm_epi_meta[,1],decreasing = T),]
write.csv(sumglm_epi_meta, file="MeanReactivitiesComp_epi_meta.csv")
x=rownames(sumglm_epi_meta)
x=unlist(stri_extract_all(x[grep(":",x)], regex = "(?<=Seq_)\\w+"))
calls_meta=pepinfo[pepinfo$peptide %in% x,]

save(vb10agL, file="var")

glm_vb10_epi_GBM=glm(Val~Seq_*Tumor, family = "gaussian", data=vb10agL, subset = D_!="meta")
sumglm_epi_GBM=summary(glm_vb10_epi_GBM)
pp=p.adjust(sumglm_epi_GBM$coefficients[,4]) 
sumglm_epi_GBM=sumglm_epi_GBM$coefficients[pp<0.05,]
sumglm_epi_GBM=sumglm_epi_GBM[order(sumglm_epi_GBM[,1],decreasing = T),]
write.csv(sumglm_epi_GBM, file="MeanReactivitiesComp_epi_GBM.csv")
x=rownames(sumglm_epi_GBM)
x=unlist(stri_extract_all(x[grep(":",x)], regex = "(?<=Seq_)\\w+"))
calls_GBM=pepinfo[pepinfo$peptide %in% x,]

ij=grep("TumorTu:", rownames(sglm))
Agj=stri_extract_all(rownames(sglm)[ij], regex="(?<=TumorTu:Ag_).+")
protosee=rownames(pepinfotab)[pepinfotab>40]
protosee=protosee[protosee %in% Agj]
nms=c("EBV BZLF","NYESO1","Claudin","HTLV","EBV env","H1_2","HLA A36","Mamaglobin","Myc","PSA","Tyrase")
names(nms)=protosee
for (p in protosee) {
  ptbl=vb10agL[vb10agL$Ag_==p,]
  glm_by_Ag=glm(Val~Seq_*D_, family = "gaussian", data=ptbl)
  sumbyAg=summary(glm_by_Ag)
  pp=p.adjust(sumbyAg$coefficients[,4]) 
  sumbyAg=sumbyAg$coefficients[pp<0.05,]
  sumbyAg=sumbyAg[order(sumbyAg[,1],decreasing = T),]
  write.csv(sumbyAg, file=paste(nms[p],".csv", sep="", collapse=""))
}

covb10=cor(t(vb10))

bmp(file="covb10.bmp", width = 10,height = 10, units = "in", res = 600)
corrplot(covb10, method = "color", order = "hclust", hclust.method = "ward.D2")
dev.off()
rm(covb10)
gc()


hclvb10=hclust(as.dist(covb10), method = "ward.D2")
pdf(file="hclvb10.pdf", width = 60, height=5)
plot(hclvb10,cex=0.1)
dev.off()

x=cutree(hclvb10, h=0.85)
table(x)

proct=proc.time()
bmp(file="covb10cut.bmp", width = 12,height = 12, units = "in", res = 600)
heatmap3(covb10[hclvb10$order,hclvb10$order], col=cpl(100), showColDendro = F,showRowDendro = F,cexRow=0.1, cexCol=0.1)
dev.off()
print(proc.time()-proct)

pcavb10=prcomp(t(vb10))
MWtnG=list()
for (i in 1:20){
    MWtnG[[i]]=wilcox.test(pcavb10$x[,i]~dgnG)
}

fviz_pca_ind(pcavb10, col.ind=dgnf, axes = 11:12)

MWtnM=list()
for (i in 1:20){
  MWtnM[[i]]=wilcox.test(pcavb10$x[,i]~dgnM)
}

fviz_pca_ind(pcavb10, col.ind=dgnf, axes = c(3,10))

MWtnC=list()
for (i in 1:20){
  MWtnC[[i]]=wilcox.test(pcavb10$x[,i]~dgnC)
}

fviz_pca_ind(pcavb10, col.ind=dgnf, axes = c(7,12))
fviz_pca_ind(pcavb10, col.ind=dgnf, axes = c(3,12))
fviz_pca_ind(pcavb10, col.ind=dgnf, axes = c(7,12))

pcalls=pcavb10$rotation[,c(3,7,10,11,12)]
hclpcalls=hclust(dist(pcalls),method="ward.D2")
heatmap3(pcalls,Rowv = as.dendrogram(hclpcalls))
fviz_pca_var(pcavb10, axes = c(7,12))
fviz_contrib(pcavb10, "var", axis=3, top = 500)
pcavb10var=get_pca_var(pcavb10)
varcor=pcavb10var$cor[,c(3,11,12)]

varcalls=apply(varcor, 2, function(x) list(x[abs(x)>quantile(abs(x),0.99)]))
x=varcalls
names(x)=NULL
vcls=unique(names(unlist(x)))
plot(cmdscale(dist(t(vb10[vcls,]))), pch=16, col=dgnf)


# IgM Reactivity Profiles of the lowest (TYR) and the highest (NY-ESO1) --------
#    biding self proteins and 4 others

bindingProfile(n=59, filename="tyrosinase") 
bindingProfile(n=3, filename="NYESO1") 
bindingProfile(n=47, filename="erbB-2")
bindingProfile(n=15, filename="HTLV1")
bindingProfile(n=38, filename="MYC")

# Search for relevant profiles using z scores -----------------------------

N=nrow(vb10)
proct=proc.time()
cl <- makeCluster(6)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cl, ex)
clusterEvalQ(cl, library(matrixStats))
clusterExport(cl, list("vb10","dgnf","N"), envir=environment())

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
    j=cbind(1:N,max.col(abs(pm)))
    pm[j]
  })
  rowSums(abs(y)>2)},  cl=cl)                                      
stopCluster(cl)                                  
print(proc.time()-proct)                         

rownames(zsets3)=rownames(vb10)
zsetmean=rowMeans(zsets3)

zsets3calls=apply(zsets3,2,function(col) rownames(zsets3)[col>1])

save(vb10,dgnG, dgnf, dgnM, dgnC, zsets3calls, file="vars")


# Recursive feature elimination --------------------------------------------
# The recursive feature elimination is performed on a HPC using the bash scripts
# cybash and runTAA as well as the R scripts runLoo.R and Loo.R in parallel
# on 21 cores. The results are passed in the data files F_1 and F_2.


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
  x=names(FSzsetst)[FSzsetst>i]
  x[order(FSzsetst[x])]
})

# Visualization in 2D using multidimensional scaling

for (i in 1:21) plotMDS(vb10[FSzscom[[i]],], col=dgnf, main=i)

# Finding level of commonality of the features that optimizes the clustring 

FSzscomcri=sapply(1:19, function(i){  # Selection of the bootstrap set 
  clucri(vb10[FSzscom[[i]],], dgnf)
})

pdf(file="hmzSetsall_99.pdf")
hm=heatmap.2(log10(vb10[FSzscom[[8]],]), hclustfun = hclwrd,  na.rm=F,key.title = NA,symkey=FALSE, cexCol=1,cexRow=0.4, colsep=c(4,10,15),rowsep=c(15,22,33,41,49,61,80,90), trace="none",col=cpl1(1000),margins = c(5,12), lwid=c(0.5,1), lhei = c(0.2,1))   #, colsep=c(4,10,15),rowsep=c(14,25,41,48,58,67,80,90)
dev.off()
pdf(file="hmzSetsall_214.pdf")
heatmap.2(vb10[FSzscom[[2]],], hclustfun = hclwrd,  na.rm=F,key.title = NA,symkey=FALSE, cexCol=1,cexRow=0.4, trace="none",col=cpl(1000), margins = c(5,12), lwid=c(0.5,1), lhei = c(0.2,1))   #, colsep=c(4,10,15),rowsep=c(14,25,41,48,58,67,80,90)
dev.off()  #  
pdf(file="hmzSetsall_329.pdf")
heatmap.2(vb10[FSzscom[[1]],], hclustfun = hclwrd,  na.rm=F,key.title = NA,symkey=FALSE, cexCol=1,cexRow=0.4, trace="none",col=cpl(1000), margins = c(5,12), lwid=c(0.5,1), lhei = c(0.2,1))   #, colsep=c(4,10,15),rowsep=c(14,25,41,48,58,67,80,90)
dev.off()
pdf(file="hmzSetsall_176.pdf")
heatmap.2(vb10[FSzscom[[3]],], hclustfun = hclwrd,  na.rm=F,key.title = NA,symkey=FALSE, cexCol=1,cexRow=0.4, trace="none",col=cpl(1000), margins = c(5,12), lwid=c(0.5,1), lhei = c(0.2,1))   #, colsep=c(4,10,15),rowsep=c(14,25,41,48,58,67,80,90)
dev.off()

# Analysis of peptide clusters --------------------------------------------
ij=cut(seq_along(FSzscom[[8]]),breaks = c(0,15,22,33,41,49,61,80,90,100), labels = F)
ctspep=aggregate(rev(FSzscom[[8]][hm$rowInd]), by=list(ij), "list")[,2]
cts=unlist(ctspep)

pats=colnames(vb10)

clstprots=lapply(ctspep,function(i) pepiused[i])

x=table(unlist(clstprots))
y=table(pepinfo$protein)
xy=cbind(y[names(x)],x)
xy[,1]=xy[,1]-xy[,2]
xn=names(x)
chsqxy=chisq.test(xy, correct = T, simulate.p.value = T)
chsqxy$stdres
1-pnorm(chsqxy$stdres[,2])
p.adjust(1-pnorm(chsqxy$stdres[,2]), method = "fdr")
xn0=xn[!xn %in% names(nms)]
nms0=c("TRP2","EBV BZLF1","EBV Capsid","CEACAM5","p53","claudin6", "HPV E7","HTLV1 env","EBNA1","G2/cyclin-B1","H1_2","HLA A36","HPV L1","HPV L2", "HHV8 LANA","Mammaglobin1","PMEL","Myc","MOG","HTLV1 gag","PSA","SSX2","erbB-2","MMP11","SOX-2","Tyrosinase","HBV X")
names(nms0)=xn0

pdf(file="callsproteinrank.pdf", width=4, height=6)
par0=par()
par(mai=c(1,1.5,0.2,0.2))
x=sort(chsqxy$stdres[,2])
barplot(x,las=2, names=nms0[names(x)],cex.names=0.8,horiz=T, xlab="Z")
par=par0
dev.off()

protscalls=lapply(ctspep,function(p) {
  x=nms0[pepinfo$protein[pepinfo$peptide %in% p]]
  names(x)=NULL
  return(x)
})

protabcls=sapply(protscalls, function(cl){
  sapply(nms0,function(n){
    sum(cl==n)
  })
})
rownames(protabcls)=nms0
viral=c(1,4,5,13,16,18,19,20,28)
colSums(protabcls[viral,])
vnms=names(nms0[viral])

callsFs=vb10[FSzscom[[8]],]
pepclsinfo=pepinfo[pepinfo$peptide %in% FSzscom[[8]],]
x=pepclsinfo$protein
names(x)=pepclsinfo$peptide
pepclsinfo=x

pdf(file="MDScallsFs.pdf", width = 5, height=5.5)
  plot(cmdscale(dist(t(callsFs))), pch=16, col=dgnf, xlab="D1", ylab="D2")
  legend("bottomright", legend=c("Control","GBM","Meta"), fill = 1:3, border=F, bty="n", cex=0.5)
dev.off()

ij=rep("self",98)
ij[pepclsinfo %in% vnms]="Vir"
calls=calls[names(pepclsinfo),]
id=rep(dgnf,each=98)
ij=rep(ij,21)
clust=sapply(names(pepclsinfo), function(p) grep(p, ctspep, fixed=T))
Seq=rep(rownames(calls), 21)
calls=data.frame(Seq=Seq,Val=c(calls),Vir=factor(ij),Dia=factor(id), Cluster=factor(clust))
pdf(file="Vir_by_cluster.pdf", width=4, height=15)
boxplot(data=calls, Val~Vir*Cluster*Dia, notch=T, las=2, horizontal=T)
boxplot(data=calls, Val~Dia*Cluster, notch=T, besides=T, horizontal=T)
dev.off()
summary(lm(data=calls, Val~(Vir+Dia+Cluster)^2))

model=aov(data=calls, Val~(Vir+Dia+Cluster)^2)
TukHSD=TukeyHSD(model, conf.level=.95)

options(max.print = 10000)
sink(file="TukeyHSDCLusters.txt")
print(TukHSD)
sink()
options(max.print = 1000)

TukHSD=lapply(TukHSD[-1], function(tb){
  tb=tb[tb[,4]<0.05,]
  tb=tb[!is.na(rownames(tb)),]
})
options(max.print = 10000)
sink(file="TukeyHSDCLusters_short.txt")
print(TukHSD)
sink()
options(max.print = 1000)

xt=sapply(xn, function(prot){
  sapply(clstprots,function(cl){
    sum(cl==prot)
  })
})
chsqxt=chisq.test(xt, simulate.p.value = T)

pnorm(t(chsqxt$stdres)) # no preferential recognition of individual proteins among clusters detectable

# GBM clusters

table(unlist(clstprots[c(1,3)]))

# Tu loss clusters
table(unlist(clstprots[c(2,5,6)]))

# GBM loss cluster

table(unlist(clstprots[4]))


# Idiotypic relations -----------------------------------------------------

load("IgJT")

XJ=IgJtrim[grep("X",IgJtrim)]
IJ=IgJtrim[-grep("X",IgJtrim)]
rm(IgJtrim)
xJs=unlist(strsplit(XJ, split="X"))
xJs=xJs[nchar(xJs)>6]
IJ=c(IJ,xJs)
nJs=nchar(IJ)
nqJ=sapply(4:14,function(n){
      sum(nJs[nJs>=n]-n+1)
})
paaIJ=colSums(alphabetFrequency(AAStringSet(IJ)))
paaIJ=paaIJ[1:20]/sum(paaIJ)

load("prots")
names(prots)=NULL
annoprots=lapply(prots, function(p) attributes(p)$Annot)
prots=sapply(prots, function(p) {
  attributes(p)=NULL
  return(p)
})
aaprots=alphabetFrequency(AAString(paste(prots, sep="", collapse="")), as.prob = T)[1:20]

nprt=nchar(prots)
nprtqs=sapply(4:14,function(n){
  sum(nprt[nprt>=n]-n+1)
})

smpls15=lapply(1:1000, function(i) sample(prots15,99))
smplsIJ=lapply(1:100, function(i) sample(IgJtrim,99))

proct=proc.time()
x=stringdistmatrix(IgJtrim,cts, "qgram")
print(proc.time()-proct)

proct=proc.time()
x=sapply(1:ncol(x),function(i) {
  table(x[,i])
})
print(proc.time()-proct)

co=1
for (p in x) {
  plot(as.numeric(names(p)),(p+0.5), ty="l", col=co, xlim=c(0,55),ylim=c(1.5,150000), log="y")
  par(new=T)
  co=co+1
}
par(new=F)


# More - 4-grams

# repertoire of 4-grams in the J regions

load("prots")
prots=unlist(prots)
q4prots=t(qgrams(prots,q=4))
q4prots=q4prots[order(q4prots, decreasing = T),]
q4IgJ=t(qgrams(IgJtrim,q=4))
q4IgJ=q4IgJ[-grep("X",rownames(q4IgJ)),]

# Bootstrapping the 4-grams distribution in J regions based on the residue 
# frequencies in the J regions - scrambled sequences

cl=makeCluster(4)
clusterExport(cl, list("IgJtrim","repScrmbl"), envir = environment())
clusterEvalQ(cl, require(stringdist))

q4BS0=pblapply(1:100, function(i){
  L=repScrmbl(IgJtrim)
  q=t(qgrams(L, q=4))
  return(q)
}, cl=cl)

stopCluster(cl)

x=c()
for (i in 1:100) x=rbind(x,q4BS0[[i]])
q4BS0=aggregate(x, by=list(rownames(x)), function(y) c(y,rep(0,100-length(y))))
nms=q4BS0$Group.1
q4BS0=q4BS0[,-1]
rownames(q4BS0)=nms
i=grep("X", nms)
q4BS0=q4BS0[-i,]
q4BStat0=cbind(rowMeans(q4BS0), rowSds(q4BS0))

q4comp0=(q4IgJ-q4BStat0[names(q4IgJ),1])/q4BStat0[names(q4IgJ),2]

q4prob0=pnorm(q4comp0)
q4sig0=q4comp0[p.adjust(q4prob0)>0.95&q4comp0>0]
q4sig0pp=unique(unlist(strsplit(names(q4sig0),split="\\.")))

# Bootstrapping the 4-grams distribution in the proteome based on the residue 
# frequencies in the J regions - scrambled sequences

cl=makeCluster(2)
clusterExport(cl, list("prots","repScrmbl"), envir = environment())
clusterEvalQ(cl, require(stringdist))

q4BSpr=pblapply(1:100, function(i){
  L=paste(sample(unlist(strsplit(prots, split=""))), sep="", collapse="")
  q=t(qgrams(L, q=4))
  return(q)
}, cl=cl)

stopCluster(cl)

x=c()
for (i in 1:100) x=rbind(x,q4BSpr[[i]])
q4BSpr=aggregate(x, by=list(rownames(x)), function(y) c(y,rep(0,100-length(y))))
nms=q4BSpr$Group.1
q4BSpr=q4BSpr[,-1]
rownames(q4BSpr)=nms
q4BStatpr=cbind(rowMeans(q4BSpr), rowSds(q4BSpr))

q4comppr=(q4prots-q4BStatpr[names(q4prots),1])/q4BStatpr[names(q4prots),2]

q4probpr=pnorm(q4comppr)
q4sigpr=q4comppr[p.adjust(q4probpr)>0.95&q4comppr>0]
q4sigprpp=unique(unlist(strsplit(names(q4sigpr),split="\\.")))

p8p=FSzscom[[8]]

cl=makeCluster(4)
clusterExport(cl, list("IgJtrim","pp","repScrmbl"), envir = environment())
clusterEvalQ(cl, require(stringdist))

q4BS=pblapply(1:100, function(i){
  L=repScrmbl(IgJtrim)
  q=t(qgrams(L, p8p, q=4))
  j=(q[,1]>0)&(q[,2]>0)
  x=q[j,1]
  names(x)=rownames(q[j,])
  return(x)
}, cl=cl)

stopCluster(cl)

q4BS=aggregate(unlist(q4BS), by=list(names(unlist(q4BS))), function(x) c(x,rep(0,100-length(x))))
nms=q4BS$Group.1
q4BS=q4BS[,-1]
rownames(q4BS)=nms
q4BStat=cbind(rowMeans(q4BS), rowSds(q4BS))

q4p8p=qgrams(IgJtrim,pp,q=4)
q4p8p=t(q4p8p)
j=(q4p8p[,1]>0)&(q4p8p[,2]>0)
x=q4p8p[j,1]
names(x)=rownames(q4p8p[j,])
q4p8p=x

q4comp=sapply(names(q4p8p), function(n){
  (q4p8p[n]-q4BStat[n,1])/q4BStat[n,2]
})

hist(q4comp, breaks=100)
q4prob=pnorm(abs(q4comp))
q4sig=q4comp[p.adjust(q4prob)>0.95&q4comp>0]
q4sigpp=unique(unlist(strsplit(names(q4sig),split="\\.")))

plot(density(log10(q4prots)), xlim = c(1,4),ylim=c(0,1.2), main="")
par(new=T)
plot(density(log10(q4prots[q4sigpp])), xlim = c(1,4), col=rgb(1,0,0,0.5), xlab="",ylim=c(0,1.2), main="")
par(new=T)
plot(density(log10(q4prots[names(q4p8p)])), xlim = c(1,4), col=rgb(0,1,0,0.5), xlab="",ylim=c(0,1.2), main="")
par(new=T)
plot(density(log10(q4prots[q4sig0pp[q4sig0pp %in% names(q4prots)]])), xlim = c(1,4), col=rgb(0,0,1,0.5), xlab="",ylim=c(0,1.2), main="")

q4names=union(names(q4IgJ),names(q4prots))
q4both=matrix(0,nrow=length(q4names), 2)
rownames(q4both)=q4names
colnames(q4both)=c("Proteins","J regions")
q4both[names(q4prots),1]=q4prots
q4both[names(q4IgJ),2]=q4IgJ
q4bth=q4both[rowProds(q4both)>0,]

q4bothsig=q4both

for (i in 1:2) q4bothsig[,i]=0
q4bothsig[q4sigprpp,1]=1
q4bothsig[q4sig0pp ,2]=1
tq4=table(q4bothsig[,1],q4bothsig[,2])
x=chisq.test(tq4)
x$stdres
tq4[2,2]/((as.numeric(sum(tq4[2,]))*as.numeric(sum(tq4[,2])))/sum(tq4))

plot(log10(q4both[q4sig0pp,]), pch=16, cex=0.5, col=rgb(0,0,0,0.5))
# J regions are represented by qgrams which correlate with a subset of the qgrams of the proteome
# the significanlty overrepresented J qgrams are better represented also in the proteome (stats?)
# the common qgrams for the calls and the J are even more overexpressed in the proteome
# the overrepresented J qgrams which are found in the calls are very highly represented in the proteome
# so - the idiotopes reflect the self antigen landscape and are the antibodies recognizing them
# change in cancer in a specific way providing diagnostic profiles.

pepsigJq4map=t(qgrams(q4sig0pp, pep, q=4))
pepsigJq4map=pepsigJq4map[rowProds(pepsigJq4map)>0,]
pepsigJq4map=rownames(pepsigJq4map)
pepsigJq4map=sapply(pep, function(p){
  q=qgrams(p,pepsigJq4map, q=4)
  q=q[,colProds(q)>0]
  print(which(pep==p))
  if (!is.null(dim(q))) res=sum(q[1,]) else res=q[[1]]
  return(res)
}) 

pepsigpq4map=t(qgrams(q4sigprpp, pep, q=4))
pepsigpq4map=pepsigpq4map[rowProds(pepsigpq4map)>0,]
pepsigpq4map=rownames(pepsigpq4map)
pepsigpq4map=sapply(pep, function(p){
  q=qgrams(p,pepsigpq4map, q=4)
  q=q[,colProds(q)>0]
  #print(which(pep==p))
  if (!is.null(dim(q))) res=sum(q[1,]) else res=q[[1]]
  return(res)
}) 

#############



########## by cluster comparisons

ctsq4sig=lapply(ctspep, function(pp) {
  x=sapply(q4sigpp, function(p) pp[grep(p, pp)])
  table(unlist(x[lengths(x)>0])) #
})

lapply(ctsq4sig, mean)
mean(unlist(ctsq4sig))

calls=cbind(calls, Idq4sig=unlist(ctsq4sig)[calls$Seq])
calls$Idq4sig[is.na(calls$Idq4sig)]=0
boxplot(calls$Idq4sig~calls$Cluster, notch=T)
boxplot(calls$Val[calls$Dia==1]~(calls$Idq4sig[calls$Dia==1]<6), notch=T)

x=aov(data=calls, Idq4sig~Cluster)
TukeyHSD(x)
# order of magnitude of No q4sig decreasing 
# 7,5,8,3,2,9,4,6,1
# *5>(3),1,8,4,6,9
# *7>3,1,8,4,6,9
# *2>3,1,8,4,6,9

x=aov(data=calls, Val~(Cluster+as.factor(Idq4sig>mean(Idq4sig))+Dia)^2)
x=TukeyHSD(x)
x=lapply(x[-2], function(y) {
  y=y[y[,4]<0.05,]
  y[order(abs(y[,1])),]  
})

q4pepall=qgrams(IgJtrim,pep,q=4)
q4pepall=t(q4pepall)
q4pepall=q4pepall[rowProds(q4pepall)>0,]

########################## to be redone

IdVirLm=lm(data=calls, Id~(Vir+Cluster)^2, notch=T)
summary(IdVirLm) # the viral epitopes to which there is 
# a uniform over- or underexpression are more idiotypically connected #
# than the respective self
boxplot(data=calls, Id~(Vir+Cluster)^2, notch=T)
IdVirAov=aov(data=calls, Id~(Vir+Cluster)^2)
IdVirAov=TukeyHSD(IdVirAov)
IdVirAov=lapply(IdVirAov[-1], function(l) {
  l=l[!is.na(l[,1]),]
  l[l[,4]<0.05,]
})
sink(file="IdVirAov.txt")
print(IdVirAov)
sink()


# Big pep comparison table ------------------------------------------------

x=quantile(q4IgJ,c(0.33,0.67))
q4pep=qgrams(IgJtrim,pep,q=4)
q4pep=t(q4pep)
j=(q4pep[,1]>0)&(q4pep[,2]>0)
q4pep=q4pep[j,1]


names(x)=names(q4pep)
q4pep=x
q4pep


# Alternative Selection Algorithm -----------------------------------------



