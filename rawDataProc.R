# "Pipeline"
# filepr preprocesses .gpr files usin gprLocNorSVM to locally 
#   normalize them by SVR, the fitted bachkgrounf is recorded 
#   as B635 columns.
# cn holds the controls' sequences
# Use pepStat to get pSet and pnSet. No spots are flagged -
#   the relative difference of the duplicates does not excede 15% in v0.
#
# v0[N spots x cases] is data from pSet, 
# v[N spots x cases] - from pnSet
# vd[N spots x (cases + flag for even chip layout rows] is 
#   reordered so that the even chip rows are in the first half 
#   and odd - in the second. This separates the duplicates in 
#   two tables of the same layout, stacked on top of each other.
# v2[N1 (peptides-controls) x 2*cases] - duplicates are now separate 
#   variables and conrtrols are removed.
# vn[sames v2] - v2 after global normalization with cyclic loess 
#   from limma, method "affy".
# vn1[N2 x cases] - means of duplicates of vn
#


require(pepStat)
require(limma)
require(Biobase)
require(sva)
require(d3heatmap)
filepr("IgMreactivities/")  # gpe folder 
mapFile="mapll.csv"
dirToParse="IgMreactivities_the"
cn=c("YPYDVPDYAG", "DYKDDDDKAS")
pSet=makePeptideSet(path=dirToParse, mapping.file = mapFile)
v0=pSet@assayData$exprs
pepts=pSet@featureRange@elementMetadata@listData$peptide
Rc=as.double(unlist(pSet@featureRange@elementMetadata@listData$row))
iRc=which(Rc %% 2==0)
pepts1=pepts[iRc]
pnSet=normalizeArray(pSet, centered = FALSE)
v=pnSet@assayData$exprs
vncol=ncol(v)
rownames(v)=pepts
rws=as.double(unlist(pSet@featureRange@elementMetadata@listData$row))%%2+1
vd=cbind(v,rws)
vd=vd[order(vd[,vncol+1]),]
v2=cbind(vd[vd[,vncol+1]==1,1:vncol],vd[vd[,vncol+1]==2,1:vncol])
v2=v2[-which(row.names(v2) %in% cn),]
pepn=rownames(v2)

pepdups=pepn[duplicated(pepn)]
i=1:length(pepdups)
pepdups2=lapply(i,function(x){which(pepn==pepdups[x])})
pepdups2=pepdups2[lengths(pepdups2)>=2]
indxy=which(duplicated(pepn))        # the same as the second column of pepdups2
v2[indxy,]=matrix(unlist(lapply(lapply(pepdups2, function(x){v2[x,]}),function(x){apply(x,2,mean)})),nrow=length(pepdups), byrow=T)
indxn=unique(c(sapply(pepdups2, function(x){x[[1]][1:(length(x[[1]])-1)]})))
indxn=indxn[!is.na(indxn)]
v2s=v2[-indxn,]
pep=rownames(v2s)

v1s=(v2s[,1:vncol]+v2s[,(vncol+1):(2*vncol)])/2
#i=c(2*(1:8)-1,17,18)   # just for the current set / 
#v=v[,i]                # remove when using a set without pre files or select only post files
vn_1=normalizeCyclicLoess(v1s[,1:10], iterations=2, method="affy")
vn_2=normalizeCyclicLoess(v1s[,11:21], iterations=2, method="affy")
v1=cbind(vn_1,vn_2)
batch=c(rep(1,10),rep(2,11))
vb1=ComBat(v1,batch, mod=dsn1, par.prior = FALSE)
vn1=normalizeCyclicLoess(vb1, iterations=1, method="affy") #vn1=vb1 

mapll=read.csv("mapll.csv",header = TRUE)
cnm=apply(as.data.frame(mapll),1,function(x){paste(c(x[3], x[5]),sep="", collapse = "_")})
colnames(vn1)=cnm
save(vn1, file="vn1")
