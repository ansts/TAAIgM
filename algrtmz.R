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
require(multcomp)
require(Biostrings)
library(lsmeans)
require(Peptides)



# Data Cleaning -----------------------------------------------------------

cpl=colorRampPalette(c("#000000FF","#0000B0FF","#00B0B0FF","#00FF0AFF","#FFFF00FF","#FF0000FF"), alpha=T)
cpl1=colorRampPalette(c("#0000AFFF","#FF0000FF"), alpha=T)
cpl2=colorRampPalette(c("#00000050","#AFFF00FF","#FF8000FF","#FF0000FF"), alpha=T)
cpl3=colorRampPalette(c("#0000FFFF","#007070FF","#0A0A0A0A","#907000FF","#FF0000FF"), alpha=T)

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
plot(rowMeans(vb10),rowSds(vb10), cex=0.3)
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
colnames(pepinfo)=c("protein","peptide")
pepinfotab=table(pepinfo[,1])
pepi=names(pepinfotab)
pepinf=pepinfo[,1]
names(pepinf)=pepinfo[,2]
pepinf=pepinf[rownames(vb10)]
pepiused=pepinf[!is.na(pepinf)]

# Correct inconsistencies in protein names

x=pepiused
x[x==pepi[[5]]]=pepi[[6]]
x[x %in% pepi[c(14,56)]]=pepi[[15]]
x[x==pepi[[24]]]=pepi[[6]]
x[x %in% pepi[c(18,19,21,22,23)]]=pepi[[25]]
x[x==pepi[[20]]]=pepi[[26]]
x[x==pepi[[35]]]=pepi[[34]]
x[x==pepi[[55]]]=pepi[[13]]
x[x==pepi[[39]]]=pepi[[41]]
x[x==pepi[[73]]]=pepi[[58]]
x[x==pepi[[44]]]=pepi[[38]]
x[x==pepi[[46]]]=pepi[[45]]
x[x==pepi[[71]]]=pepi[[68]]

xinfotab=table(x)
pepi=names(xinfotab)
pepiused=x
viral=c(2,4,9:18,20,25:32,37:39,42,45:47,49,52)

# General Table for ANOVA -------------------------------------------------

load("shortNames")
vb10ag=vb10[names(pepiused),]-mean(vb10[names(pepiused),])
vb10agL=melt(vb10ag)
colnames(vb10agL)=c("Seq_","D_","Val")
x=stri_extract_all(vb10agL$D, regex="(?<=_)\\w+" )
vb10agL$D_=unlist(x)
x=pepiused[vb10agL$Seq_]
vb10agL=cbind(vb10agL,Ag_=x)
vb10agL$D_=as.factor(vb10agL$D_)
vb10agL$Ag_=as.factor(vb10agL$Ag_)
x=aggregate(vb10agL$Val, by=list(vb10agL$Ag_), "mean")
x[which.min(abs(scale(x[,2]))),1]

Y=unique(vb10agL$Ag_)
Aglm=t(sapply(Y, function(a){
  x=vb10agL[vb10agL$Ag_==a,]
  x=lm(data = x, Val~D_)
  x=summary(x)
  fv=x$fstatistic
  p=pf(fv[1],fv[2],fv[3],lower.tail=F)
  xeff=x$coefficients[,1]
  names(xeff)=paste("Eff",names(xeff), sep="_")
  x=x$coefficients[,4]
  names(x)=paste("p.val",names(x), sep="_")
  return(c(xeff[2:3],x[2:3],Overall.p=p))
}))
rownames(Aglm)=Y
Aglm[,5]=p.adjust(Aglm[,5],method="BH")
sink(file="SigAgReact.txt")
print(Aglm[Aglm[,5]<0.1,])
sink()

TuSqG=TukeySq(vb10agL, "GBM") 
x=TuSqG$coefficients[,1]
p=TuSqG$coefficients[,4]

table(Coeff=x>0,p=p<0.05, pepiused)

TuSqC=TukeySq(vb10agL, "cntr") 
x=TuSqC$coefficients[,1]
p=TuSqC$coefficients[,4]

table(Coeff=x>0,p=p<0.05, pepiused)

TuSqM=TukeySq(vb10agL, "meta") 
x=TuSqM$coefficients[,1]
p=TuSqM$coefficients[,4]

table(Coeff=x>0,p=p<0.05, pepiused)

vb10agL$Ag_=shortNames[vb10agL$Ag_]

Daov=aov(data=vb10agL, Val~D_)
TukeyHSD(Daov, "D_")

# IgM Reactivity Profiles ------- 
#    of the lowest (TYR) and the highest (NY-ESO-1) biding self proteins
# as well as some other

bindingProfile(n=59, filename="tyrosinase") 
bindingProfile(n=3, filename="NYESO1") 
bindingProfile(n=50, filename="erbB-2")
bindingProfile(n=15, filename="HTLV1")
bindingProfile(n=49, filename="endoHTLV1")
HLA36ag=bindingProfile(n=24, filename="HLA_A36")

