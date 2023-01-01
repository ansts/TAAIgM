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
require(nlme)


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
pats=colnames(vb10)
plot(rowMeans(vb10),rowSds(vb10), cex=0.3)
pdf(file="vb10dotplots.pdf", width=10, height=10)
rng=range(vb10)
plot(as.data.frame(vb10),cex=0.1, pch=16, xlim=rng, ylim=rng)
dev.off()
save(vb10, file="vb10")
save(list=c("pats","dgn","dgnf","dgnC","dgnG","dgnM"), file="ColnamesAndDiag")

# blood group
bgr=c(ez="0",gj="0",idmtr="0",id="A",ls="B",zs="A",gr="0",yt="0",rp="B",
      lb="0",vs="A",ri="A",md="A",ai="0",dk="A",yk="B",he="0",nn="0",ab="AB",hs="0",jd="0")
patnm=colnames(v0)
patnm=unlist(stri_extract_all(patnm, regex="\\w+(?=_)"))
bgr=bgr[patnm]
BlGr=data.frame(A=grepl("A",bgr), B=grepl("B",bgr), N=grepl("0",bgr))
rownames(BlGr)=names(bgr)
save(list=c("sex","BlGr"), file="BloodGroupAndSex")


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
save(list=c("pep","pepinfo","pepi","pepiused"), file="pepbioinfo")

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
meanAg=aggregate(rowMeans(vb10ag), by=list(pepiused), "mean")
meanAg=meanAg[order(meanAg$x),]
vb10agL=melt(vb10ag)
colnames(vb10agL)=c("Seq_","Pat_","Val")
x=stri_extract_all(vb10agL$Pat_, regex="(?<=_)\\w+" )
vb10agL=cbind(vb10agL,D_=unlist(x))
x=pepiused[vb10agL$Seq_]
vb10agL=cbind(vb10agL,Ag_=x)
vb10agL$D_=as.factor(vb10agL$D_)
vb10agL$Pat_=as.factor(vb10agL$Pat_)
vb10agL$Ag_=factor(vb10agL$Ag_, levels=meanAg$Group.1, ordered = T)

Aglm=lme(Val~D_*Ag_, random=~1|Pat_, data = vb10agL)
res=lsmeans(Aglm, pairwise~D_|Ag_)
reres=summary(res$contrasts)
j=reres$p.value<0.1
rres=data.frame(Contrast=as.character(reres$contrast[j]),Antigen=as.character(reres$Ag_[j]),Estimate=as.numeric(reres$estimate[j]),p.value=as.numeric(reres$p.value[j]))
rres=rres[order(rres$Estimate),]
sink(file="SigAgReact.txt")
  print(reres)
sink()

resAg=summary(lsmeans(Aglm,~D_*Ag_, random=~1|Pat_))

pdf(file="resAg.pdf", width=10, height=50)
  plot(resAg)
dev.off()


# IgM Reactivity Profiles ------- 
#    of the lowest (TYR) and the highest (NY-ESO-1) biding self proteins
#    as well as some other

bindingProfile(n=59, filename="tyrosinase") 
bindingProfile(n=3, filename="NYESO1") 
bindingProfile(n=50, filename="erbB-2")
bindingProfile(n=15, filename="HTLV1")
bindingProfile(n=49, filename="endoHTLV1")
HLA36ag=bindingNNProfile(n=24, filename="HLA_A36")

