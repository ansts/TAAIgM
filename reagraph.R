require(igraph)
require(reshape2)
require(parallel)
require(pbapply)
require(stringdist)
require(matrixStats)
require(chisq.posthoc.test)
require(ppcor)
require(stringi)
require(dunn.test)
require(limma)
require(gplots)
require(nlme)
require(multcomp)

require(irr)
require(pROC)
require(uwot)
require(RColorBrewer)
require(qualV)
require(infotheo)
require(ggplot2)
require(rgl)
require(corrplot)
require(Peptides)
require(Biostrings)
require(factoextra)
require(FactoMineR)
require(riverplot)
require(class)

# require(statnet)

load("IgJtrim")
viral=c(2,4,9:18,20,21,26:35,41,42,45,46,48:50,52,55,60)
vir=pepinfo$peptide[pepinfo$protein %in% pepi[viral]]

# Adjacency Matrix --------------------------------------------------------

vs10=scale(vb10, scale = F)
vb010=vb10-min(vb10)+1

ij=t(combn(length(pep),2))
  
J=cut(seq_along(ij[,1]), 1000, labels = F)
cl=makeCluster(4)
clusterExport(cl, c("ij","J","vs10","valF"))

Fs=pbsapply(1:max(J), function(i){
      ii=which(J==i)   
      v=apply(ij[ii,],1,function(j){
        x=valF(vs10[j,])
        return(x)
      })
      return(v)
}, cl=cl)
stopCluster(cl)
Fs=unlist(Fs, recursive = F)

save(vs10,pep,file="vars")

ij=rbind(ij,cbind(ij[,2],ij[,1]))
x=data.frame(X=ij[,1], Y=ij[,2], V=c(Fs,Fs))
ajmFs=acast(x,X~Y)
rownames(ajmFs)=pep
colnames(ajmFs)=pep

pp=pepinfo[pepinfo$peptide %in% rownames(ajmFs),]
x=pp$protein
names(x)=pp$peptide
rm(pp)

adjmdiag=t(sapply(seq_along(pep)[-1], function(j) {
  d=stringdist(pep[j],pep[j-1], method = "lcs")
  print(c(pep[j],pep[j-1],d,ajmFs[j,j-1]))
  return(c(d,ajmFs[j,j-1])) 
}))

adjmsmpl=sapply(1:10000, function(i){
  fl=T
  while(fl){
    ij=sample(seq_along(pep),2)
    d=stringdist(pep[ij[1]],pep[ij[2]], method="lcs")
    fl=abs(diff(ij))<12|d<25
  }  
  ajmFs[ij[1],ij[2]]
})
hist(log10(adjmsmpl), xlim=c(-1.5,2), breaks = 50)
par(new=T)
hist(log10(adjmdiag[adjmdiag[,1]<7,2]), col=rgb(1,0,0,0.5),xlim=c(-1.5,2), breaks = 50, border = F)

d1=density(log10(adjmsmpl), adjust=0.5)
d2=density(log10(adjmdiag[adjmdiag[,1]<7,2]), adjust=0.5)
plot(d1, lwd=2,ylim=c(0,1.5), xlim=c(-1.5,1.5), main="", xlab="log10 F value")
par(new=T)
plot(d2, lwd=2,ylim=c(0,1.5), xlim=c(-1.5,1.5),col=2,  xlab="", ylab="", main="")
par(new=T)
plot(c(log10(3),log10(3)),c(0,1), ty="l",col=8,ylim=c(0,1.5), xlim=c(-1.5,1.5), xlab="",ylab="")

df=data.frame(Group=c(rep(1,length(adjmsmpl)), rep(2,length(adjmdiag[adjmdiag[,1]<7,2]))),Reac=c(log10(adjmsmpl),log10(adjmdiag[adjmdiag[,1]<7,2])))
reaRoc=roc(data=df, Group, Reac, auc=T)
plot(reaRoc,print.auc=T, print.thres="best", ci=T)
plot(reaRoc,print.auc=T, print.thres=log10(3), ci=T)

#  Graph ------------------------------------------------------------------

Gvs=graph_from_adjacency_matrix(ajmFs*(ajmFs>3), weighted = T, diag = F, mode = "undirected")
Gvs=simplify(Gvs)
components(Gvs)

ccGvs=transitivity(Gvs, type=("localundirected"))
ccGvs[is.na(ccGvs)]=0



names(Jq4map)=names(V(Gvs))
Jq4map[names(pepsigJq4map)]=pepsigJq4map

pq4map=rep(0,length(names(V(Gvs))))
names(pq4map)=names(V(Gvs))
pq4map[names(pepsigpq4map)]=pepsigpq4map

pp=pepiused[names(pepiused) %in% vGvs]
linkseqs=vGvs[!(vGvs %in% names(pepiused))]
x=union(union(grep("SGSGS", linkseqs, value=T),grep("GSGSGS", linkseqs, value=T)),grep("GSGSG", linkseqs, value=T))
linker=rep("linker",length(x))
names(linker)=x
y=setdiff(linkseqs,x)
marg=rep("marg",length(y))
names(marg)=y
pp=c(pp,linker,marg)
pp=pp[vGvs]
prtsn=as.factor(pp)
prts=unlist(as.numeric(prtsn))
prtlbl=levels(as.factor(pp))

mx=vb10-min(vb10)+1

Gvlets=graphlets(Gvs, rep(1,length(E(Gvs))))

bgr=c(ez="0",gj="0",idmtr="0",id="A",ls="B",zs="A",gr="0",yt="0",rp="B",
      lb="0",vs="A",ri="A",md="A",ai="0",dk="A",yk="B",he="0",nn="0",ab="AB",hs="0",jd="0")
patnm=colnames(v0)
patnm=unlist(stri_extract_all(patnm, regex="\\w+(?=_)"))
bgr=bgr[patnm]
BlGr=data.frame(A=grepl("A",bgr), B=grepl("B",bgr), N=grepl("0",bgr))
rownames(BlGr)=names(bgr)
chisq.posthoc.test(table(sex,BlGr$A), simulate.p.value=T)
chisq.posthoc.test(table(sex,BlGr$B), simulate.p.value=T)
Abgr=BlGr$A+BlGr$N
Bbgr=BlGr$B+BlGr$N
chisq.posthoc.test(table(sex,Bbgr), simulate.p.value=T)

# Graph pars -------------------------------------------------------------

Gvs=set_vertex_attr(Gvs, "protein", value = prts)
Gvs=set_vertex_attr(Gvs, "calls", value = 1*(vGvs %in% FSzscom[[8]]))
Gvs=set_vertex_attr(Gvs, "IgJsig", value =Jq4map )
Gvs=set_vertex_attr(Gvs, "Protsig", value =pq4map[names(V(Gvs))])
Gvs=set_vertex_attr(Gvs, "GvsC", value =GvsC[names(V(Gvs))])
Gvs=set_vertex_attr(Gvs, "MvsC", value =MvsC[names(V(Gvs))])
Gvs=set_vertex_attr(Gvs, "Mean", value =rowMeans(vs10)[names(V(Gvs))])
Gvs=set_vertex_attr(Gvs, "CV", value =rowSds(vb010[names(V(Gvs)),])/rowMeans(vb010)[names(V(Gvs))])

virCluster=names(ego(Gvs,1,"DHILEPSIPWKSK")[[1]])
for (p in virCluster){
  plot(vs10[p,order(dgnf)], ty="l", ylim=c(-3,4))
  par(new=T)
}
par(new=F)

secCluster=names(ego(Gvs,1,"GTLVALVGLFVLLAF")[[1]])
for (p in secCluster){
  plot(vs10[p,order(dgnf)], ty="l", ylim=c(-3,4))
  par(new=T)
}
par(new=F)


backCluster=names(ego(Gvs,1,"GLCKAIQEQCCFLN")[[1]])
for (p in backCluster){
  plot(vs10[p,order(dgnf)], ty="l", ylim=c(-3,4))
  par(new=T)
}
par(new=F)

back2Cluster=names(ego(Gvs,1,"NDVCAQVHPQKVTKF")[[1]])
for (p in back2Cluster){
  plot(vs10[p,order(dgnf)], ty="l", ylim=c(-3,4))
  par(new=T)
}
par(new=F)

forCliquer=names(ego(Gvs,1,"VPVSEPVPEPEPEPE")[[1]])
for (p in forCliquer){
  plot(vs10[p,order(dgnf)], ty="l", ylim=c(-3,4))
  par(new=T)
}
par(new=F)


# NN Stats of Diagnosis ---------------------------------------------------
cl=makeCluster(15)
clusterExport(cl, c("Gvs","vs10","dgnf"))
clusterEvalQ(cl, require(igraph))
NNdiaov=t(pbsapply(pep,function(p) {
  pp=names(ego(Gvs,1, p)[[1]])
  x=vs10[pp,]
  if (!is.null(nrow(x)))  {
    n=nrow(x)
    Y=c(x)
    pat=as.factor(rep(1:21,each=n))
    xA=as.factor(rep(dgnf,each=n))
    iG=xA!=3
    iM=xA!=2
    pAG=summary(aov(Y[iG]~xA[iG]+Error(pat[iG])))[[1]][[1]][[5]][1]
    pAM=summary(aov(Y[iM]~xA[iM]+Error(pat[iM])))[[1]][[1]][[5]][1]
    dAG=diff(aggregate(Y[iG], by=list(xA[iG]), "mean")$x)
    dAM=diff(aggregate(Y[iM], by=list(xA[iM]), "mean")$x)
    xA=c(`2-1`=dAG,`3-1`=dAM,pG=pAG,pM=pAM)
    return(xA)
  }
  else {
    pAG=wilcox.test(x[dgnf!=3]~dgnf[dgnf!=3])[[3]][[1]]
    pAM=wilcox.test(x[dgnf!=2]~dgnf[dgnf!=2])[[3]][[1]]
    dAG=diff(aggregate(x[dgnf!=3], by=list(dgnf[dgnf!=3]), "mean")$x)
    dAM=diff(aggregate(x[dgnf!=2], by=list(dgnf[dgnf!=2]), "mean")$x)
    xA=c(`2-1`=dAG,`3-1`=dAM,pG=pAG,pM=pAM)
    return(xA)
  }
},cl=cl))
stopCluster(cl)

rownames(NNdiaov)=pep
NNdiaov[,3:4]=apply(NNdiaov[,3:4],2,p.adjust,method="BH",n=3*length(NNdiaov[,4]))

for (j in 1:3) {
  NNdiaov[(NNdiaov[,j+3]>0.001),j]=0
}

Gvs=set_vertex_attr(Gvs, "GvsCkw", value =NNdiaov[,1])
Gvs=set_vertex_attr(Gvs, "MvsCkw", value =NNdiaov[,2])
write.graph(Gvs, format = "graphml", file="Gvs.graphml")

mx=vb10-min(vb10)+1
NNdimn=t(sapply(pep,function(p) {
  pp=names(ego(Gvs,1, p)[[1]])
  x=mx[pp,]
  print(which(vGvs==p))
  if (!is.null(nrow(x))) n=nrow(x) else n=1
  Y=c(x)
  pat=as.factor(rep(1:21,each=n))
  xA=aggregate(aggregate(Y,by=list(pat), "mean")[,2], by=list(dgnf), "mean")[,2]
  return(c(log2(xA[2]/xA[1]),log2(xA[3]/xA[1]),log2(xA[3]/xA[2])))
}))
colnames(NNdimn)=c("2-1", "3-1","3-2")

NNdicmb=NNdimn*(abs(NNdimn)>0.1)*(NNdiaov[,4:6]<0.001)
Gvs=set_vertex_attr(Gvs, "Gvscmb", value =NNdicmb[,1])
Gvs=set_vertex_attr(Gvs, "Mvscmb", value =NNdicmb[,2])
write.graph(Gvs, format = "graphml", file="Gvs.graphml")

table(cut(NNdicmb[,1], c(-1e6,-1e-6,1e-6,1e6), labels=c("Lost","None","Gained")))
table(cut(NNdicmb[,2], c(-1e6,-1e-6,1e-6,1e6), labels=c("Lost","None","Gained")))
table(cut(NNdicmb[,3], c(-1e6,-1e-6,1e-6,1e6), labels=c("Meta","None","GBM")))
sum(NNdicmb[,1]!=0)/nrow(NNdicmb)
sum(NNdicmb[,2]!=0)/nrow(NNdicmb)
sum(NNdicmb[,3]!=0)/nrow(NNdicmb)


Zmod(Gvs, 1+(NNdicmb[,1]!=0))
Zmod(Gvs, 1+(NNdicmb[,2]!=0))
Zmod(Gvs, 1+(NNdicmb[,3]!=0))
Zmod(Gvs, 1+(NNdicmb[,1]>0))
Zmod(Gvs, 1+(NNdicmb[,2]>0))
Zmod(Gvs, 1+(NNdicmb[,3]>0))
Zmod(Gvs, 1+(NNdicmb[,1]<0))
Zmod(Gvs, 1+(NNdicmb[,2]<0))
Zmod(Gvs, 1+(NNdicmb[,3]<0))

Zass(Gvs, 1+(NNdicmb[,1]!=0))
Zass(Gvs, 1+(NNdicmb[,2]!=0))
Zass(Gvs, 1+(NNdicmb[,3]!=0))
Zass(Gvs, 1+(NNdicmb[,1]>0))
Zass(Gvs, 1+(NNdicmb[,2]>0))
Zass(Gvs, 1+(NNdicmb[,3]>0))
Zass(Gvs, 1+(NNdicmb[,1]<0))
Zass(Gvs, 1+(NNdicmb[,2]<0))
Zass(Gvs, 1+(NNdicmb[,3]<0))


# Cliques -----------------------------------------------------------------


clqGvs=max_cliques(Gvs, min=3) 
tclql=table(lengths(clqGvs))
vgrep=Vectorize(grep, vectorize.args = "pattern")
ls=as.numeric(names(tclql))

mx=vb10-min(vb10)+1

ij=cut(seq_along(clqGvs), 100, labels=F)
cl=makeCluster(15)
clusterExport(cl, c("clqGvs","ij", "mx","dgnf"))

diaclqall=pblapply(1:100, function(i){
  jj=which(ij==i)
  res=t(sapply(jj, function(j){
    cq=clqGvs[[j]]
    cq=names(cq)
    y=c(mx[cq,])
    x=as.factor(rep(dgnf,each=length(cq)))
    data.frame(Stat=kruskal.test(y~x)$statistic, P=kruskal.test(y~x)$p.value)
  }))
  return(res)
}, cl=cl)
stopCluster(cl)

x=c()
for (i in 1:100){
  x=rbind(x,diaclqall[[i]])
}

diaclqall=as.data.frame(apply(x,2,unlist))

diaclqall$P=p.adjust(diaclqall$P, method = "BH")

posclqall=clqGvs[diaclqall$P<0.05]
tposclqall=table(lengths(posclqall))

clqord=clqGvs[order(lengths(clqGvs))]
x=list()
xp=c()
for (i in seq_along(clqGvs)) {
  y=names(clqord[[i]])
  n=sum(y %in% xp)
  if (n==0) x[[i]]=y else x[[i]]=c()
  xp=union(xp,y)
}
x=x[lengths(x)>0]

uniclq=x
y0=unique(unlist(uniclq))
res=t(sapply(1:1000,function(j){
  y1=sample(y0,40)
  yr1=clucri(vs10[y1,], dgnf)
  y2=sample(pep[!(pep %in% y0)],40)
  yr2=clucri(vs10[y2,], dgnf)
  return(cbind(yr1,yr2))
}))

hist(res[,1], xlim=range(res))
par(new=T)
hist(res[,2], xlim=range(res), col=2)

hist(dGvs[uniclq])
pepbycl=table(unlist(lapply(clqGvs, names)))
pepbyclpo=table(unlist(lapply(posclqall, names)))


# Cliques 2 ---------------------------------------------------------------

clqGvs10=max_cliques(Gvs, min=3, max=10)
mx=vb10-min(vb10)+1

diaclq=as.data.frame(t(sapply(clqGvs10, function(cq){
  cq=names(cq)
  y=c(mx[cq,])
  x=as.factor(rep(dgnf,each=length(cq)))
  data.frame(Stat=kruskal.test(y~x)$statistic, P=kruskal.test(y~x)$p.value)
})))

diaclq$P=p.adjust(diaclq$P, method = "BH")

posclq=clqGvs10[diaclq$P<0.1]
tposclq=table(lengths(posclq))

barplot(tclql, xlab= "Clique Size", ylab="N", yaxt="n", xaxt="n")
axis(1, cex.axis=0.75, at=(0:35)*1.2+0.7, labels=3:38, las=2)
axis(2, cex.axis=0.75)
barplot(tposclqall, xlab= "Clique Size", ylab="N", yaxt="n", xaxt="n")
axis(1, cex.axis=0.75, at=(0:35)*1.2+0.7, labels=3:38, las=2)
axis(2, cex.axis=0.75)

barplot(tposclqall/tclql[1:length(tposclqall)], 
        xlab= "Clique Size", ylab="Proprtion with DE in Diagnoses", yaxt="n", xaxt="n")
axis(1, cex.axis=0.75, at=(0:30)*1.2+0.7, labels=3:33, las=2)
axis(2, cex.axis=0.75)

clucriposcl=sapply(posclq, function(x) clucri(vs10[x,], dgnf))
boxplot(clucriposcl~lengths(posclq), notch=T)

proct=proc.time()
clucri3posclq=sapply(seq_along(posclq), function(i){
  if ((i %% 1000)==0) {
    print(i)
    print(proc.time()-proct)
  }
  x=posclq[[i]]
  x=names(x)
  ij=combn(length(x),3)
  median(apply(ij,2,function(j) clucri(vs10[x[j],], dgnf)))
})
boxplot(clucri3posclq~lengths(posclq), notch=T)

x=lapply(posclq, names)
x=sapply(x, function(y){
  cbind(seq=y,n=as.numeric(length(y)))
})
x0=c()
for (l in x){
  x0=rbind(x0,l)
}
x0=as.data.frame(x0)
x0[,2]=as.numeric(x0[,2])
x=aggregate(x0[,1], by=list(x0[,2]), function(y) list(unique(y)))
clcrposclbysize=sapply(seq_along(x$x), function(i){
  l=x$x[[i]]
  res=t(sapply(1:1000,function(j){
    y=sample(l,40)
    print(c(i,j))
    yr=clucri(vs10[y,], dgnf)
    cbind(Crit=yr, D=mean(log10(dGvs[y])), L=x$Group.1[i])
  }))
  return(res)
})



x0=c()
for (clm in 1:8){
  y=clcrposclbysize[,clm]
  y=cbind(y[1:1000],y[1001:2000], y[2001:3000])
  
  x0=rbind(x0,y)
}

clqposBSclcrXD_ClqSz=x0
z=aggregate(clqposBSclcrXD_ClqSz[,2], by=list(clqposBSclcrXD_ClqSz[,3]), "mean")[,2]
z=rep(z, each=1000)
boxplot(clqposBSclcrXD_ClqSz[,1]~clqposBSclcrXD_ClqSz[,3], notch=T, xlab="Clique Size", ylab="Clustering Criterion")
boxplot(clqposBSclcrXD_ClqSz[,2]~clqposBSclcrXD_ClqSz[,3], notch=T, xlab="Clique Size", ylab="Mean Log(Degree)")
cpl3=colorRampPalette(c("#00000080","#0000F080","#00C00980","#00FF0080","#FFFF0080","#FF000080"), alpha=T)
for(i in x$Group.1){
  plot(clqposBSclcrXD_ClqSz[clqposBSclcrXD_ClqSz[,3]==i,2],
       clqposBSclcrXD_ClqSz[clqposBSclcrXD_ClqSz[,3]==i,1], 
       col=cpl3(10)[i], pch=16, cex=0.3,
       xlim=range(clqposBSclcrXD_ClqSz[,2]),
       ylim=range(clqposBSclcrXD_ClqSz[,1]))
  #xy=lm(clqposBSclcrXD_ClqSz[clqposBSclcrXD_ClqSz[,3]==i,1]~clqposBSclcrXD_ClqSz[clqposBSclcrXD_ClqSz[,3]==i,2])
  #abline(xy, col=i)
  par(new=T)
}
par(new=F)

plot()

summary(lm(clqposBSclcrXD_ClqSz[,1]~clqposBSclcrXD_ClqSz[,2]+clqposBSclcrXD_ClqSz[,3]))

cor(x[,1],clcrposclbysize[,2])

# Diagnosis clustering correlates inversely with the log degree of the peptides.
# It does not dependent on the clique size independent of the degree.
# Clique size correlates, of course, very strongly with the mean log degree. 


callsposclq=unique(names(unlist(posclq)))
dClq=sapply(clqGvs10, function(cq){
  sum(dGvs[names(cq)])-length(cq)
})
diaclq=cbind(diaclq, D=dClq)

x1=log10(unlist(diaclq$Stat))
x2=log10(unlist(diaclq$D)+0.5)
y=data.frame(St=x1, D=x2)
summary(lm(data=y, log10(St+0.5)~D))

diaclq=cbind(diaclq, ClqSz=lengths(clqGvs10))

boxplot(data=diaclq,log10(unlist(Stat))~ClqSz)

x=dGvs[callsposclq]
x1=names(x)[x<25]
#cls=rfe(vs10[x1, ], dgnf)
cls=Loolocal0(x1,vs10)
#x2=rownames(cls[[1]])
x2=cls
x3=names(x)[!(names(x) %in% x1)]

x0=c(x3, x1[!(x1 %in% x2)])
best=c()
maxcr=0
repeat {
  cri=0
  p0=c()
  for (p in x0){
    x4=c(x2,p)
    cr=clucri(mx[x4,], dgnf)
    if (cr>cri) {
      cri=cr
      p0=p
    }
  }
  x2=c(x2,p0)
  cr0=clucri(mx[x2,], dgnf)
  if (cr0>maxcr){
    maxcr=cr0
    best=x2
  }
  print(cr0)
  x0=setdiff(x0,p0)
  if (length(x0)==0) break
}

plotMDS(vs10[best,], col=dgnf)


# Feature selection Loo cliques 2 -----------------------------------------

#callsclqfKW=callsclqf

Gi=lapply(1:21, function(i){
  l0=l[[1]]
  for (l1 in l[-1]) l0=rbind(l0,l1)
  M=array(0,dim=c(length(pp), length(pp)))
  M[l0[,1:2]]=l0[,3]
  M[l0[,c(2,1)]]=l0[,3]
  rownames(M)=pp
  colnames(M)=pp
  G=graph_from_adjacency_matrix(M, mode="undirected",weighted = T, diag=F)
  G=simplify(G)
  return(G)
})

N=nrow(vb10)
mx=vb10-min(vb10)+1 
J=cut(seq_along(clqGvs3_4), 45, labels=F)
 
clqFs=lapply(1:21, function(pat){
  print(pat)
  proct=proc.time()
  clqGvs3_4=max_cliques(Gi[[pat]], min=3, max=4)
 
  cl=makeCluster(15)
  clusterExport(cl, c("pat","mx","dgnf","N", "clucri","connfix","dunnfix","bhgamfix","BHgamma","clqGvs3_4","J"), envir = environment())
  diaclq=pbsapply(1:45, function(j){
    iJ=J==j
    sapply(clqGvs3_4[iJ], function(cq){
          cq=names(cq) 
          n=length(cq)
          y=mx[cq,-pat]
          x=clucri(y, dgnf[-pat])
          bsmx=array(sample(mx[,-pat]), dim = dim(mx[,-pat]))
          bsx=sapply(1:1000,function(i) {
            bscq=sample(1:N,n)
            clucri(bsmx[bscq,],dgnf[-pat])
          })
          f=ecdf(bsx)
          x=f(x)
          return(x)
        })
      }, cl=cl)
  stopCluster(cl)

  diaclq=unlist(diaclq)
  posclq=clqGvs3_4[diaclq>0.95]
  callsposclq=unique(names(unlist(posclq)))
  print(length(callsposclq))
  x1=callsposclq
  cls=Loolocal0(x1,vs10)
  return(cls)
})

cofmxclqFs=SVMforValidation(clqFs,X=mx,dgnf )
clqFTb=table(unlist(clqFs))
mx=vb10-min(vb10)+1
for (i in 0:20){
  y=names(clqFTb)[clqFTb>i]
  plotMDS(mx[y,], col=dgnf, main=i)
  print(c(i,clucri(vs10[y,], dgnf)))
}
y=names(clqFTb)[clqFTb>4]
plotMDS(mx[y,], col=dgnf, labels=NULL, pch=16)
callsclqf=names(clqFTb)[clqFTb>4]

pdf(file="hmclqFall_66.pdf")
hm=heatmap.2(vs10[callsclqf,], hclustfun = hclwrd,  na.rm=F,key.title = NA,symkey=FALSE, cexCol = 0.7,cexRow = 0.5,
             trace="none",col=cpl1(1000),colsep=c(6,10,15),rowsep=c(7,13,21,28,38,46,52,59),
             margins = c(5,12), lwid=c(0.5,1), lhei = c(0.2,1))    
dev.off()
clcalls=callsclqf[rev(hm$rowInd)]
patord=hm$colInd
cli=cut(seq_along(clcalls),breaks = c(0,7,13,21,28,38,46,52,59,70), labels = F)
pati=cut(1:21,breaks = c(0,6,10,15,22), labels = F)
claov=sapply(1:9,function(i){
  n=length(cli[cli==i])
  
  y=c(mx[clcalls[cli==i],patord])
  cl=as.factor(rep(pati,each=n))
  p=as.factor(rep(1:21,each=n))
  a=aov(y~cl/p)
  x=summary(a)
  ta=TukeyHSD(a,"cl")
  ta=ta$cl[,1]*(ta$cl[,4]<0.01)
  return(ta)
})
closs=apply(claov,2,function(j) any(j[c(3,5,6)]>0))
clain=apply(claov,2,function(j) any(j[c(3,5,6)]<0))

x=clcalls
x[x=="GMALRVTRNSKINAE"]="ALRVTRNSKINAENK"
protclcls=sapply(x, function(p){
  pepinfo[pepinfo[,2]==p,1]
})
tclcls=table(unlist(protclcls),cli)
chisq.test(tclcls, simulate.p.value = T)
chsqph_prot_cl=chisq.posthoc.test(tclcls, simulate.p.value = T)
x=as.matrix(chsqph_prot_cl[chsqph_prot_cl$Value=="p values",3:11])
i=rowSums(x<0.05)
chsqph_prot_cl[which(i==1)*2,]
j=colSums(x<0.05)
virclcls=clcalls %in% vir

tvirclcls=table(virclcls,cli)
chisq.posthoc.test(tvirclcls, simulate.p.value=T)
tvirloss=table(virclcls,cli %in% 4:5)
chisq.posthoc.test(tvirloss, simulate.p.value=T)
tvirgain=table(virclcls,cli %in% 6:7)
chisq.posthoc.test(tvirgain, simulate.p.value=T)
# Increased reactivity to viruses correlates with tumors - 
# EBV with GBM and papilloma vir with all

ctspep=aggregate(clcalls, by=list(cli), "list")[,2]

clij=rep(0,length(vGvs))
names(clij)=vGvs
for (i in 1:9){
  clij[ctspep[[i]]]=i
}

Gvs=set_vertex_attr(Gvs, "CallsCl", value=clij)

write.graph(Gvs, format = "graphml", file="Gvs.graphml")

clsIg7tb=table(Calls=as.factor(clij>0), Id=as.factor(q7ppp>0))
chisq.posthoc.test(clsIg7tb, simulate.p.value=T)
clqtbcut=cut(clqFTb, c(0,1.5,2.5,20), labels=F)
clqFIg7tb=table(Fs=as.factor(vGvs %in% names(clqFTb[clqFTb==1])), Id=as.factor(q7ppp>0))
chisq.posthoc.test(clqFIg7tb, simulate.p.value=T)

boxplot(log10(dGvs)~(clij>0), notch=T)
kruskal.test(log10(dGvs)~clij)

chisq.posthoc.test(table(cut(AB0mNNg[,1], c(-1e3,-1e-6,1e-6,1e3), labels = F),clij>0), simulate.p.value=T)
chisq.posthoc.test(table(cut(AB0mNNg[,2], c(-1e3,-1e-6,1e-6,1e3), labels = F),clij>0), simulate.p.value=T)
chisq.posthoc.test(table(cut(AB0mNNg[,3], c(-1e3,-1e-6,1e-6,1e3), labels = F),clij>0), simulate.p.value=T)
chisq.posthoc.test(table(cut(AB0mNNg[,3], c(-1e3,-1e-6,1e-6,1e3), labels = F),clij), simulate.p.value=T)
# The calls avoid reactivities cross reactive with BlGr A 
# The highest expression calls (cluster 1) have also higher expression in men 
# then non calls

# Mapping DEx Cliques  ----------------------------------------------------

clq345Gvs=max_cliques(Gvs, min=3, max=5)
mx=vb10-min(vb10)+1

diaclq345=as.data.frame(t(sapply(clq345Gvs, function(cq){
  cq=names(cq)
  y=c(mx[cq,])
  x=as.factor(rep(dgnf,each=length(cq)))
  data.frame(Stat=kruskal.test(y~x)$statistic, P=kruskal.test(y~x)$p.value)
})))

diaclq345$P=p.adjust(diaclq345$P, method = "holm")

posclq34=clq34Gvs[diaclq34$P<0.05]
posclq345pep=table(unlist(lapply(posclq345, names)))


diaclqDnT=t(sapply(posclq345,function(cq) {
  x=names(cq)
  y=c(mx[x,])
  x=as.factor(rep(dgnf,each=length(x)))
  d=dunn.test(y,g=x, method="bh")
  dz=d$Z
  names(dz)=d$comparisons
  return(dz)  #cbind(d$Z,d$P.adjusted)
}))
rownames(diaclqDnT)=seq_along(diaclqDnT[,1])

diaclqDnT[abs(diaclqDnT)<1.65]=0
clqDnT=sign(diaclqDnT)

peptoclq=lapply(seq_along(posclq34), function(i){
  clq=posclq34[[i]]
  data.frame(pep=names(clq), clq=i)
})
x=c()
for (tb in peptoclq){
  print(tb)
  x=rbind(x,tb)
}
peptoclq=x
peptoclq=aggregate(peptoclq$clq, by=list(peptoclq$pep), "list")
x=t(apply(peptoclq,1,function(l){
  if(length(l[2][[1]])==1) {
     x=clqDnT[unlist(l[2]),]
    } else x=colSums(clqDnT[unlist(l[2]),])
  return(x)
}))
rownames(x)=peptoclq[,1]
peptoclq=x
pclqf=sign(peptoclq)
pclqf[abs(peptoclq)<0.5]=0

for (i in 0:7){
  x=names(clqFTb)[clqFTb>i]
  z=rownames(pclqf)
  y=union(x,z)
  print(i+1)
  print(table(loo=as.factor(y %in% x), all=as.factor(y %in% z)))
}

x=unique(pclqf)
x=cbind(x,x[,1]+2*x[,2]-4*x[,3])
x[order(x[,4]),]
y=pclqf[,1]+2*pclqf[,2]-4*pclqf[,3]
i=cut(y, c(-6,-4,-2.5,-1.5,0,1.5,2.5,3.5), labels=F)
li=c("Gu","GMu","Md","Gmd","gMu","Mu","GMd")
cbind(x[order(x[,4]),],li)
y=li[i]
names(y)=rownames(pclqf)
clqflab=rep("N", length(vGvs))
names(clqflab)=vGvs
clqflab[names(y)]=y

Gvs=set_vertex_attr(Gvs, "ClqFsLab", value=clqflab)

write.graph(Gvs, format = "graphml", file="Gvs.graphml")



######### Comparisons -----------------------------------------------------
pp=pepinfo$peptide[pepinfo$peptide %in% pep]
calls=vs10[pp,]
x=pepinfo$protein
names(x)=pepinfo$peptide
x=x[pp]
pepclsinfo=x
vnms=read.table("vir.txt")[,1]
ij=rep("self",length(pp))
ij[pepclsinfo %in% vnms]="Vir"

q7pp=qgrams(IgJtrim,pp,q=7)
q7pp=t(q7pp)
q7pp=q7pp[(q7pp[,1]>0)&(q7pp[,2]>0),]
  
q7pp=stringdistmatrix(pp,rownames(q7pp), method = "lcs")
ls=abs(nchar(pp)-7)
ls=array(rep(ls,ncol(q7pp)), dim=dim(q7pp))
q7pp=rowSums(q7pp==ls)

cls=pp %in% FSzscom[[8]]
cls0=pep %in% FSzscom[[8]]
names(cls0)=pep
sigprt=pepsigpq4map[pp]

Deg=degree(Gvs)
Deg=Deg[pp]
ClC=transitivity(Gvs,type="localundirected")
names(ClC)=names(V(Gvs))
ClC=ClC[pp]
ClC[is.na(ClC)]=0

x=data.frame(calls,Pept=rownames(calls),vir=as.factor(ij),protein=as.factor(pepclsinfo), Id=q7pp, cls=as.factor(cls), SigSeq=as.factor(sigprt), Degree=Deg, ClCoeff=ClC)
calls=melt(x, measure.vars = 1:21, variable.name = "Patient")
Dia=unlist(stri_extract_all(calls$Patient, regex="(?<=_)\\w+"))
calls=cbind(calls, Dia=as.factor(Dia))

boxplot(data=calls,value~Dia+SigSeq, las=2, xlab="", notch=T)
Taovy=lapply(levels(calls$Dia), function(D){
  aovy=aov(data=calls,value~SigSeq, subset=Dia==D)
  boxplot(data=calls,value~SigSeq, subset=Dia==D, main=D)
  x=TukeyHSD(aovy, "SigSeq")
  x$SigSeq[x$SigSeq[,4]<0.05, ]
})

summary(glm(data=calls, value~(cls+Id+log10(Degree+0.5))^2, family = gaussian))

chsqClsId=chisq.test(table(as.factor(cls),q7pp), simulate.p.value = T)
chsqClsId$stdres
chsqClsSigP=chisq.test(table(as.factor(cls),sigprt>11), simulate.p.value = T)
chsqClsSigP$stdres
chsqClsDeg=chisq.test(table(as.factor(cls),Deg>22), simulate.p.value = T)
chsqClsDeg$stdres

summary(glm(cls~(log10(ClC+0.01)+log10(Deg+0.5))^2, family = binomial))
boxplot(data=calls, value~Id+cls, notch=T, las=2, xlab="")

ticks=(1:8)*2-0.5
boxplot(data=calls, value~cls+Id, col=rep(c(2,3),8), outline=F, las=2, xlab="Number of Overlapping Id 7-mers", ylab="Binding", xaxt='n')
axis(side=1,ticks,labels = c(0:6,8))
id=2*calls$Id
abline(lm(data=calls, value~id, subset = cls==TRUE), col=3, lwd=2)
abline(lm(data=calls, value~id, subset = cls==FALSE), col=2, lwd=2)
legend("bottom", c("Non-calls", "Calls"), fill=c(2,3), bty="n", cex=0.75, y.intersp = 2)

boxplot(data=calls, log10(Degree+0.5)~cls, outline=F,ylab="Log10 Degree", xlab="",names=c("Non-calls", "Calls"), notch=T)

x=cbind(x,value=rowMedians(as.matrix(x[,1:21])))
plot(x$Degree[x$cls==FALSE]+0.5,x$value[x$cls==FALSE], xlim=range(x$Degree+0.5), ylim=range(x$value), cex=0.5, log="x", pch=16, col=rgb(0.9,0.1,0,0.03), xlab="Degree", ylab="Binding")
par(new=T)
plot(x$Degree[x$cls==TRUE]+0.75,x$value[x$cls==TRUE], xlim=range(x$Degree+0.5), ylim=range(x$value), cex=0.5, log="x", pch=16,  col=rgb(0,0.5,0,0.6), xlab="Degree", ylab="Binding")
abline(lm(data=x, value~log10(Degree+0.5), subset = cls==TRUE), col=3, lwd=3)
abline(lm(data=x, value~log10(Degree+0.5), subset = cls==FALSE), col=2, lwd=3)
legend("bottom", c("Non-calls", "Calls"), fill=c(2,3), bty="n", cex=0.75, y.intersp = 2)

boxplot(data=x, log10(ClCoeff+0.5)-0.301~cls, notch=T, ylab="Log10 Clustering Coefficient", names=c("Non-calls","Calls"), xlab="")
legend("topright",legend="p=0.068", bty="n")
summary(lm(data=x, log10(ClCoeff+0.5)~cls,))

plot(log10(Deg+0.5)+rnorm(length(Deg), 0,0.01),log10(clqmmbr[pp]+0.5)+rnorm(length(pp), 0,0.01), pch=16, cex=cls*0.5+1, col=cls*1+1)
x=log10(ClC[Deg>3]+0.01)+rnorm(length(ClC[Deg>3]), 0,0.005)
y=log10(clqmmbr[pp][Deg>3]+0.5)+rnorm(length(pp[Deg>3]), 0,0.03)
plot(x,y, xlim=range(x), ylim=range(y), pch=16, cex=0.3)
par(new=T)
plot(x[cls[Deg>3]==TRUE],y[cls[Deg>3]==TRUE],col=rgb(1,0,0,1), pch=16, xlim=range(x), ylim=range(y))

xy=table(CC=x>(-0.5),Calls=cls[Deg>3]==TRUE)
chsqxy=chisq.test(xy, simulate.p.value = T)
chsqxy$stdres
CCxClqn_c=table(CC=x>(-0.3),ClqN=y>2, Calls=cls[Deg>3]==TRUE)
CCxClqn_c=as.matrix(CCxClqn_c)
CCxClqn_c=cbind(CCxClqn_c[1:4],CCxClqn_c[5:8])
chsqxy=chisq.test(CCxClqn_c, simulate.p.value = T)
chsqxy$stdres

boxplot(y~cls[Deg>3][x>-0.3], notch=T)
summary(lm(y~(x+cls[Deg>3])^2))
plot3d(x,y,log10(Deg[Deg>3]), col = as.numeric(cls[Deg>3])+1, size =5)

Grnd=erdos.renyi.game(length(V(Gvs)),length(E(Gvs)),type="gnm")
ClCrnd=transitivity(Grnd,type="localundirected")
ClCrnd[is.na(ClCrnd)]=0
Degrnd=degree(Grnd)
Clqrnd=max_cliques(Grnd, min = 4)

pdf(file="Clqmem_Deg.pdf", width = 15, height = 15)
plot(log10(Deg[Deg>3]+0.5),y, cex=0)
text(log10(Deg[Deg>3]+0.5),y, cex=0.3, labels = pp)
dev.off()

weirdclqmem=c("SPWAPKKHRRLSSDQ","YKTCKQAGTCPPDII","SNFKVRDPVIQERLD","ETTDLYCYEQLNDSS","GDICNTMHYTNWTHI","GGAGAGGAGAGGGGR","SGSSITTLSSGPHSL","VQKRHYLVFLIAHHY","RMGGIQRRDHLQTLR","PEPLPQGQLTAYHVS","YPLLKLLGSTWPTTP","LKLLGSTWPTTPPRP","GTHGGTGAGAGAGGA","PQKRPSCIGCKGTHG","APARRTRKPQQPESL")

q7pps=qgrams(IgJtrim,names(V(Gvs)),q=7)
q7pps=t(q7pps)
q7pps=q7pps[(q7pps[,1]>0)&(q7pps[,2]>0),]

q7ppp=stringdistmatrix(names(V(Gvs)),rownames(q7pps), method = "lcs", useNames = T)
ls=abs(nchar(names(V(Gvs)))-7)
ls=array(rep(ls,ncol(q7ppp)), dim=dim(q7ppp))
q7ppp=rowSums(q7ppp==ls)

q7mult=rownames(q7pps[q7pps[,1]>1,])
IgJmult=lapply(q7mult, function(p){
  grep(p,IgJtrim, value = T)
})
IgJmult=unique(unlist(IgJmult))
seqtoFASTA(IgJmult, "IgJmult.txt")

ppq7mult=unique(unlist(lapply(q7mult, function(p){
  grep(p, pep, value = T)
})))

ijppq7=stri_locate_all_fixed(ppq7mult,"GSGS")
ijppq7GS=sapply(ijppq7, function(x) x[1])
names(ijppq7GS)=ppq7mult
ijppq7GS=ijppq7GS[!is.na(ijppq7GS)]
ij=cut(ijppq7GS, c(0,3.5,7.5,16), labels=F)
for (i in 1:3){
  for (p in names(ijppq7GS)[ij==i]){
    plot(vb10[p,], ty="l", ylim=range(vb10), main=i)
    par(new=T)
  }
  par(new=F)
}


# New comp ----------------------------------------------------------------

x=sapply(4:38, function(n){
  pprev=unique(unlist(lapply(clqGvs[lengths(clqGvs) %in% 3:(n-1)], names)))
  p=unique(unlist(lapply(clqGvs[lengths(clqGvs)==n], names)))
  setdiff(p,pprev)
})
boxplot(sapply(x[1:7], function(l) dGvs[l]))

# cliques larger than 10 are formed by recombination only, 
# cliques up to 5 contain 99.2% of the peptides
# the mean degree increases to 6 cliques
lengths(x)

# More graph par ----------------------------------------------------------


diameter(Gvs, weights=rep(1, length(E(Gvs)))) #  12, with weights - 40.69
Zmod(Gvs, vertex.attributes(Gvs)$protein+1)
Zmod(Gvs, vertex.attributes(Gvs)$Ig7+1)
Zmod(Gvs, vertex.attributes(Gvs)$Vir+1)
pwc=power_centrality(Gvs, exponent=-.5)


Zass(Gvs, vertex.attributes(Gvs)$Mean)
Zass(Gvs, vertex.attributes(Gvs)$protein)
Zass(Gvs, vertex.attributes(Gvs)$Protsig)
Zass(Gvs, vertex.attributes(Gvs)$Ig7)
Zass(Gvs, vertex.attributes(Gvs)$Vir)

ClCg=transitivity(Gvs, type = "localundirected")
ClCg[is.na(ClCg)]=0
ego.size(Gvs,)

Gvs=set_vertex_attr(Gvs, "Ig7", value=q7ppp)
Gvs=set_vertex_attr(Gvs, "CC", value=ClCg)
Gvs=set_vertex_attr(Gvs, "PwC", value=pwc)

write.graph(Gvs, format = "graphml", file="Gvs.graphml")

Gid=induced.subgraph(Gvs, V(Gvs)[vertex.attributes(Gvs)$Ig7>0])
components(Gid)
# from 564 Id 7-mer hits 464 are in a single component and 82 in 68 components (only one of 4 vertices)

tDeg=table(Deg)
plot(log10(as.numeric(names(tDeg))),log10(tDeg))
ld=Deg+0.5
dj=cut(ld,seq(range(ld)[1]*0.9,range(ld)[2]+1.1,diff(range(ld))/50), labels = F)
tDegl=aggregate(ld, by=list(dj), "length")
plot(tDegl, log="xy")
plot(tDegl[,1],tDegl[,2], pch=16, cex=0.5, xlab="Degree", ylab="N", log="y")
lmtDegl=lm(log10(tDegl[,2])~tDegl[,1])
abline(lmtDegl)
text(37,316, labels="y=473*0.884^x")
text(37,230, labels="R^2=0.956, p<1e-15")

# Louvain clustering ------------------------------------------------------

vGvs=names(V(Gvs))
GvLouv=treeCrawLouv(vGvs, Gvs) #recursive Louvain clustering - yields a tree
ni=sapply(GvLouv, function(y) length(unique(lengths(y))))
i=which(ni>1)
x=lapply(i,function(j) lapply(GvLouv[[j]], unlist))
x=unlist(x, recursive = F)
GvLouv=c(GvLouv[-i], x)
x=GvLouv[lengths(GvLouv)>4]
x=lapply(x, unlist)
GvLouv=x  # Final clustering of the whole graph
GvcL=GvLouv
names(GvcL)=seq_along(GvcL)
GvcL=lapply(GvcL, function(l) {
  l=unlist(l)
  names(l)=NULL
  return(l)
})
GvcL=melt(GvcL)
x=GvcL$L1
names(x)=GvcL$value
GvcL=x

#  Map to the graph -------------------------------------------------------

Gvs=set_vertex_attr(Gvs, "Louv", value=GvcL[vGvs])
write.graph(Gvs, format = "graphml", file="Gvs.graphml")


# Graph of clusters -------------------------------------------------------

cx=components(Gvs)
GVscon=induced.subgraph(Gvs, names(cx$membership[cx$membership==1]))

clmx=sapply(1:14, function(i){
      print(i)
      sapply(1:14, function(j){
        if (j==i) {
          return(0)
        }
        else {
          l=Lvcl[Lvcl %in% c(i,j)]
          v=names(l)
          g=induced_subgraph(GVscon, v)
          modularity(g, l[names(V(g))]+1)
        }
      })
})
xy=cmdscale(clmx)
rownames(xy)=0:14
plot(xy, cex=0)
text(xy, labels=0:14)
hist(clmx)
clmxin=1/(clmx+0.5)
diag(clmxin)=0
GLv=graph_from_adjacency_matrix(clmxin, mode="undirected", weighted = T)
write.graph(GLv, format = "graphml", file="GLv.graphml")
plot3d(embed_laplacian_matrix(GLv, no=3)$X, size = 5)

plot3d(embed_adjacency_matrix(GVscon, no=10)$X)
# PCA shows 271 dimensional space for 90% of the variation

vGVscon=names(V(GVscon))
Lvcntr=Lvcl[vGVscon]
GLvcntr=contract(GVscon,mapping=Lvcntr, vertex.attr.comb = "ignore")
GLvcntr=simplify(GLvcntr, edge.attr.comb = "sum")
write.graph(GLvcntr, format = "graphml", file="GLvc.graphml")

# Louvain clusters correlations -------------------------------------------

# Correlation with calls

tGcsLC=table(GvcL,names(GvcL) %in% FSzscom[[8]])

chisqGvsLC=chisq.test(tGcsLC, simulate.p.value = T)
chisqGvsLC$stdres
pcls=p.adjust(pnorm(abs(chisqGvsLC$stdres[,2]), lower.tail = F))

# Correlation with proteins 

pp=pepinfo$protein
names(pp)=pepinfo$peptide
tGcsLCp=table(GvcL, pp[names(GvcL)])
chisqGvsLCpp=chisq.test(tGcsLCp+0.5, simulate.p.value = T)
ppp=array(p.adjust(pnorm(abs(chisqGvsLCpp$stdres), lower.tail = F)), dim=dim(chisqGvsLCpp$stdres))
ij=which(ppp<0.05, arr.ind = T)
resLC=cbind(rownames(chisqGvsLCpp$stdres)[ij[order(ij[,1]),1]],colnames(chisqGvsLCpp$stdres)[ij[order(ij[,1]),2]], chisqGvsLCpp$stdres[ij[order(ij[,1]),]])


x=names(GvcL)[GvcL==5][names(GvcL)[GvcL==5] %in% FSzscom[[8]]]
z=sapply(x,function(y) {
  j=grep(y,clstprots)
  print(j)
  clstprots[[j]][names(clstprots[[j]][!is.na(clstprots[[j]])])==y]
})

Cl5=names(GvcL)[col==5]

LuvClLCS8=lapply(0:14, function(i){
  print(i)
  Cl=names(GvcL)[col==i]
  ClLCS=sapply(Cl, function(p1){
    p1=unlist(strsplit(p1, split=""))
    sapply(Cl, function(p2){
      p2=unlist(strsplit(p2, split=""))  
      l=LCS(p1,p2)$LCS
      paste(l, sep="", collapse="")
    })
  })
  ClLCS=table(ClLCS[lower.tri(ClLCS)])
  ClLCS=cbind(nchar(names(ClLCS)),ClLCS)
  print(ClLCS)
  ClLCS8=names(ClLCS[ClLCS[,1]==8,2])
  return(ClLCS8)
})
seqtoFASTA(LuvClLCS8[[5]], "Clust4r5LCS.txt")

Cl8=names(col)[col==8]

Cl5H12=Cl5[Cl5 %in% pepinfo$peptide[pepinfo$protein==pepi[[23]]]]
Cl5BZ=Cl5[Cl5 %in% pepinfo$peptide[pepinfo$protein==pepi[[2]]]]
Cl8Tyr=Cl8[Cl8 %in% pepinfo$peptide[pepinfo$protein==pepi[[59]]]]

strmx=as.matrix(stringdistmatrix(c(Cl5H12,Cl5BZ,Cl8Tyr), method = "lcs"))
image(strmx<14)

# Correlation with idiotopes

q7cut=cut(q7ppp, c(-1,0.5,1.5,100), labels = F)
ti7GLC=table(col,q7cut)

chisqi7GLC=chisq.posthoc.test(ti7GLC, simulate.p.value = T)

pi7=p.adjust(pnorm(abs(chisqi7GLC$stdres), lower.tail = F))
ij=which(array(pi7,dim=dim(chisqi7GLC$stdres))<0.05, arr.ind = T)
ij=cbind(rownames(chisqi7GLC$stdres)[ij[,1]],colnames(chisqi7GLC$stdres)[ij[,2]])
ij=ij[order(ij[,1]),]

bincalli7sLC=cbind(calls=tGcsLC[,2]/rowSums(tGcsLC), I7=(rowSums(ti7GLC)-ti7GLC[,1])/rowSums(ti7GLC))
barplot(t(bincalli7sLC), beside = T)
barplot(t(cbind(tGcsLC[,1],tGcsLC[,2])), xlab = "Clusters", ylab="No of sequences", las=2)
legend(13.4,700, legend=c("Calls", "Non-calls"), cex=0.85, bty = "n", fill=c("Light Gray", "Black"))
text(6.5,700,"*", cex=2)
barplot(t(cbind(ti7GLC[,1],ti7GLC[,2],ti7GLC[,3])), xlab = "Clusters", ylab="No of sequences", las=2)
legend(1,740, legend=c(">1 IgJ Sequences", "  1 IgJ Sequences", "  0 IgJ Sequences"),cex=0.85, bty = "n", fill=c("Light Gray", "Dark Gray","Black"))
text(c(5.5,6.7,10.3),c(420,360,744),"*", cex=1.4)
barplot(rowSums(tGcsLCp>0), ylim = c(0,60), las=2)
barplot(t(tGcsLCp))

rs=rowSums(tGcsLCp)/sum(tGcsLCp)
names(rs)=NULL
cs=colSums(tGcsLCp)/sum(tGcsLCp)
names(cs)=NULL
x0=rs %*% t(cs)
x=tGcsLCp/sum(tGcsLCp)
xl=x
xl[xl>0]=log(xl[xl>0])
x0l=x0
x0l[x0l>0]=log(x0l[x0l>0])
mi=x*(xl-x0l)
image(mi, col=cpl1(100))

# Correlation with viral

tvGLC=table(GvcL,vGvs %in% vir)

chisqvGLC=chisq.posthoc.test(tvGLC, simulate.p.value = T)

pepvircl6=vGvs[GvcL==6 & vGvs %in% vir]
length(unique(pepinfo$protein[pepinfo$peptide %in% pepvircl6]))
length(vnms)


# Correlations with NN Dia Groups -----------------------------------------

for (gr in colnames(NNdiaov[,1:3])){
  y=NNdiaov[names(Lvcl),gr]
  cy=cut(y, c(-1e6,-1e-6,1e-6,1e6), labels=c("Negative", "Zerro","Positive"))
  print(table(cy,Lvcl))
  print(chisq.posthoc.test(table(cy,Lvcl), simulate.p.value=T))
}

tbNNLv=table(rowSums(NNdiaov[,1:2])!=0, Lvcl)
csqtbNNLv=chisq.posthoc.test(tbNNLv, simulate.p.value=T)
csqtbNNLv=as.matrix(apply(csqtbNNLv[c(1,3),3:17],2,as.numeric))
barplot(csqtbNNLv, beside=T ,col=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), xlab="Cluster", ylab="Standardized Residuals")
lines(c(0,46),c(3,3))
lines(c(0,46),c(-3,-3))
legend("topleft", legend=c("Under represented","Over represented"), fill=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), bty="n")

tbNNGLv=table(NNdiaov[NNdiaov[,1]!=0,1]>0, Lvcl[NNdiaov[,1]!=0])
csqtbNNGLv=chisq.posthoc.test(tbNNGLv, simulate.p.value=T)
csqtbNNGLv=as.matrix(apply(csqtbNNGLv[c(1,3),3:ncol(csqtbNNGLv)],2,as.numeric))
barplot(csqtbNNGLv, beside=T ,col=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), xlab="Cluster", ylab="Standardized Residuals", las=2, main="Lost in GBM")
lines(c(0,46),c(3,3))
lines(c(0,46),c(-3,-3))
legend("topleft", legend=c("Under represented","Over represented"), fill=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), bty="n")
tbNNGLv=table(NNdiaov[NNdiaov[,1]!=0,1]<0, Lvcl[NNdiaov[,1]!=0])
csqtbNNGLv=chisq.posthoc.test(tbNNGLv, simulate.p.value=T)
csqtbNNGLv=as.matrix(apply(csqtbNNGLv[c(1,3),3:ncol(csqtbNNGLv)],2,as.numeric))
barplot(csqtbNNGLv, beside=T ,col=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), xlab="Cluster", ylab="Standardized Residuals", las=2, main="Over Expressed in GBM")
lines(c(0,46),c(3,3))
lines(c(0,46),c(-3,-3))
legend("topleft", legend=c("Under represented","Over represented"), fill=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), bty="n")
tbNNMLv=table(NNdiaov[NNdiaov[,2]!=0,2]>0, Lvcl[NNdiaov[,2]!=0])
csqtbNNMLv=chisq.posthoc.test(tbNNMLv, simulate.p.value=T)
csqtbNNMLv=as.matrix(apply(csqtbNNMLv[c(1,3),3:ncol(csqtbNNMLv)],2,as.numeric))
barplot(csqtbNNMLv, beside=T ,col=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), xlab="Cluster", ylab="Standardized Residuals", las=2, main="Lost in Meta")
lines(c(0,46),c(3,3))
lines(c(0,46),c(-3,-3))
legend("topleft", legend=c("Under represented","Over represented"), fill=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), bty="n")
tbNNMLv=table(NNdiaov[NNdiaov[,2]!=0,2]<0, Lvcl[NNdiaov[,2]!=0])
csqtbNNMLv=chisq.posthoc.test(tbNNMLv, simulate.p.value=T)
csqtbNNMLv=as.matrix(apply(csqtbNNMLv[c(1,3),3:ncol(csqtbNNMLv)],2,as.numeric))
barplot(csqtbNNMLv, beside=T ,col=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), xlab="Cluster", ylab="Standardized Residuals", las=2, main="Over Expressed in Meta")
lines(c(0,46),c(3,3))
lines(c(0,46),c(-3,-3))
legend("topleft", legend=c("Under represented","Over represented"), fill=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), bty="n")

HLA36=names(pp)[pp==pepi[24]]

HLA36=NNdiaov[HLA36,1:2]
HLA36pos=HLA36[rowSums(HLA36)!=0,]
HLA36pep=rownames(HLA36pos)

cyclB1=names(pp)[pp=="G2/mitotic-specific cyclin-B1 (UniProt P14635)"]

cyclB1=NNdiaov[cyclB1,1:2]
cyclB1pos=cyclB1[rowSums(cyclB1)!=0,]
cyclB1pep=rownames(cyclB1pos)

PSA=names(pp)[pp==pepi[40]]

PSA=NNdiaov[PSA,1:2]
PSApos=PSA[rowSums(PSA)!=0,]
PSApep=rownames(PSApos)

Myc=names(pp)[pp==pepi[35]]

Myc=NNdiaov[Myc,1:2]
Mycpos=Myc[rowSums(Myc)!=0,]
Mycpep=rownames(Mycpos)

p53nn=names(pp)[pp=="Cellular tumor antigen p53 (UniProt Q2XN98)"]

p53nn=NNdiaov[p53nn,1:2]
p53nnpos=p53nn[rowSums(p53nn)!=0,]
p53nnpep=rownames(p53nnpos)

H4=names(pp)[pp=="Histone H4 (UniProt P62805)"]

H4nn=NNdiaov[H4,1:2]
H4nnpos=H4nn[rowSums(H4nn)!=0,]
H4nnpep=rownames(H4nnpos)

Tyrnn=names(pp)[pp=="Tyrosinase (UniProt P14679)"]

Tyrnn=NNdiaov[Tyrnn,1:2]
Tyrnnpos=Tyrnn[rowSums(Tyrnn)!=0,]
Tyrnnpep=rownames(Tyrnnpos)

x=c(HLA36pep, PSApep, Mycpep, cyclB1pep, p53nnpep, H4nnpep, Tyrnnpep)
plotMDS(vs10[x,], col=dgnf, pch=16, cex=3)
clucri(vb[x,], dgnf)

l=x
L=list(dgnC, dgnM, dgnG)
cl=makeCluster(3)
clusterExport(cl, c("rfe", "BHgamma", "bhgamfix","connfix","dunnfix","vb10","L","dgnf","l"), envir=environment())
calls=pblapply(L, function(d) {
  m=vb10[l,]
  x=rfe(m, d)
  rownames(x[[1]])
},cl=cl)
x=table(unlist(calls))
fs=names(x)[x>1]
stopCluster(cl)

plotMDS(vb10[fs,], pc=16, col=dgnf, cex = 3)


# HLA A 36

HLA36=HLA36[-1:2,]
plot(HLA36[,2], ty="l")
plot(HLA36)
tbHLA36=table(GBM=cut(HLA36[,1], c(-1e6,-1e-6,1e-6,1e6), 
                  labels=c("Negative", "Zerro","Positive")),
              Meta=cut(HLA36[,2], c(-1e6,-1e-6,1e-6,1e6), 
                  labels=c("Negative", "Zerro","Positive")))
chisq.posthoc.test(tbHLA36, simulate.p.value=T)

# Correlations with sex ---------------------------------------------------


sex=c("M","M","F","F","F","F","M","F","M","F","M","F","F","M","M","M","M","M","F","F","F")
names(sex)=pats
Lvcl=col
names(Lvcl)=rownames(vb10)

wilsex=sapply(1:14, function(i){
  print(i)
  m=vb10[Lvcl==i,]
  m=melt(m)
  m=cbind(m, Sex=sex[m$Var2])
  w=wilcox.test(data=m, value~Sex)
  return(c(W=w$statistic, p=w$p.value))
})
wilsex[2,]=p.adjust(wilsex[2,])

ij=c(5,6,9,10,12)
ry=range(vb10[rownames(vb10) %in% names(Lvcl)[Lvcl %in% c(5,6,9,10,12)],])
boxplot(vb10[rownames(vb10) %in% names(Lvcl)[Lvcl %in% ij]]~Lvcl[Lvcl %in% ij])
for (i in ij) {
  hist(vb10[rownames(vb10) %in% names(Lvcl)[Lvcl==i],], breaks=50, xlim=ry, main=i)
}

# Cor. Sex - NN - egonetwork ----------------------------------------------


y=sapply(vGvs,  function(n){
  x=names(ego(Gvs, nodes=n)[[1]])
  if (length(x)<2) return(NULL)
  x=melt(vb10[x,])
  x=cbind(x,sex=sex[x$Var2])
  wilcox.test(data=x,value~sex)$p.value
})
y=unlist(y[!is.na(y)])
y=p.adjust(y)
sum(y<0.05)
snn=vGvs %in% names(y)[y<0.05]
names(snn)=vGvs
sxnn=sapply(vGvs, function(n){
  x=names(ego(Gvs, nodes=n)[[1]])
  if (length(x)<2) return(NULL)
  x=melt(vb010[x,])
  x=cbind(x,sex=sex[x$Var2])
  rt=aggregate(value~sex, data=x, FUN="median")
  rt=rt$value[rt$sex=="F"]/rt$value[rt$sex=="M"]
  log10(summary(aov(value~sex, data=x))[[1]][1,4])*sign(log10(rt))
})
sxnn=unlist(sxnn)
sxn=rep(0,length(vGvs))
names(sxn)=vGvs
sxn[names(sxnn)]=sxnn
sxn[names(snn)[!snn]]=0

table(M=sxn<1, F=sxn>-1)

sxn_nn=rep(0, length(vGvs))
names(sxn_nn)=vGvs
jf=unique(unlist(sapply(names(sxn)[sxn>0], function(v){
  names(ego(Gvs,1, nodes=v)[[1]])
})))
jm=unique(unlist(sapply(names(sxn)[sxn<0], function(v){
  names(ego(Gvs,1, nodes=v)[[1]])
})))
jfm=intersect(jf,jm)
jf=setdiff(jf,jfm)
jm=setdiff(jm,jfm)

sxn_nn[jf]=1
sxn_nn[jm]=-1
sxn_nn[jfm]=2
Gvs=set_vertex_attr(Gvs, "Sexnbh", value=sxn)
Gvs=set_vertex_attr(Gvs, "SexNN", value=sxn_nn)
write.graph(Gvs, format = "graphml", file="Gvs.graphml")

tbclsxnn=table(cls0,sxn_nn)
chsqtbclsxnn=chisq.posthoc.test(tbclsxnn, simulate.p.value = )
chsqtbclsxnn

jfmset=read.csv(file="JFM.csv")
sort(table(pepinfo$protein[pepinfo$peptide %in% jfmset[,1]]), decreasing = T)
jfmvset=read.csv(file="JFMvir.csv")
sort(table(pepinfo$protein[pepinfo$peptide %in% jfmvset$Label]), decreasing = T)

tsxvir=table(sxn_nn, vGvs %in% vir)
chisq.posthoc.test(tsxvir, simulate.p.value=T)
jj=(sxn>0)+1+(sxn<0)*2
jj[jj==1]="none"
jj[jj==2]="F"
jj[jj==3]="M"
tsxvir1=table(vGvs %in% vir,jj)
chisq.posthoc.test(tsxvir1, simulate.p.value=T)
# Gender specific reactivities correlate with antiviral in women
# and are anticorrelated in men

tsxi7=table(sxn_nn,q7cut)
chsqtsxi7=chisq.posthoc.test(tsxi7, simulate.p.value=T) 
barplot(as.matrix(chsqtsxi7[chsqtsxi7$Value =="Residuals",3:5]), beside=T, 
        col=c(rgb(0,0,1,1),rgb(0.5,0.5,0.5,0.5),rgb(1,0,0,1), rgb(1,1,0,1)),
        names=c("no Id", " 1 Id 7-frame", ">1 Id 7-frame"), 
        ylab="Chi Square Residuals", ylim=c(-5,5))
text(c(6.5,8.5,14.5),rep(4,3),labels = c("*","*","*"), cex=3)

# Highly common idiotope determinants are crossreactive with both male and 
# female specific reactivties (yellow) but are not indifferent to gender 
# separataion. Unique idiotype reactivities are more common in man than in
# women.

tsxi7_=table((sxn>-1)-(sxn<1),q7cut)
chsqtsxi7_=chisq.posthoc.test(tsxi7_, simulate.p.value=T) 

# When using only the central vertex instead of the whole neighborhood
# these findings do not reach statistical significance

sxn_n=sxn_nn
sxn_n[sxn_n==2]=0
tsx1i7=table(sxn_n,q7cut)
chsqtsx1i7=chisq.posthoc.test(tsx1i7, simulate.p.value=T) 
barplot(as.matrix(chsqtsx1i7[chsqtsx1i7$Value =="Residuals",3:5]), beside=T, 
        col=c(rgb(0,0,1,1),rgb(0.5,0.5,0.5,0.5),rgb(1,0,0,1)),
        names=c("no Id", " 1 Id 7-frame", ">1 Id 7-frame"), 
        ylab="Chi Square Residuals", ylim=c(-5,5))
text(c(5.5,7.5),rep(4,2),labels = c("*","*"), cex=1.5)
text(c(11.5,11.5),c(4,3.3),labels = c("\U2640","\U2642"), col=c(rgb(1,0,0,1),rgb(0,0,1,1)),cex=2)


tviri7=table(vGvs %in% vir, q7cut)
legend("bottomright", 
       fill=c(rgb(0,0,1,1),rgb(0.5,0.5,0.5,0.5),rgb(1,0,0,1), rgb(1,1,0,1)),
       legend = c("M","none","F","FM"))

chisq.posthoc.test(tviri7, simulate.p.value=T) 

tnnsx=table(rowSums(NNdiaov[,1:2])!=0,sxn_nn)
chsqtnnsxn=chisq.posthoc.test(tnnsx, simulate.p.value=T)
chsqtnnsxn=as.matrix(apply(chsqtnnsxn[c(1,3),3:6],2,as.numeric))
barplot(chsqtnnsxn, beside=T ,col=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), 
        xlab="Gender Related Compartments", ylab="Standardized Residuals",
        names=c("M", "None", "F","Overlap"))
lines(c(0,46),c(3,3))
lines(c(0,46),c(-3,-3))
lines(c(0,46),c(2.5,2.5), col="grey")
lines(c(0,46),c(-2.5,-2.5), col="grey")
legend("topleft", legend=c("Deficient","Abundant"), fill=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), bty="n")

tnnGsx=table(NNdiaov[NNdiaov[,1]!=0,1]<0,sxn_nn[NNdiaov[,1]!=0])
chsqtnnGsxn=chisq.posthoc.test(tnnGsx, simulate.p.value=T)
chsqtnnGsxn=as.matrix(apply(chsqtnnGsxn[c(1,3),3:6],2,as.numeric))
barplot(chsqtnnGsxn, beside=T ,col=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), 
        xlab="Gender Related Compartments", ylab="Standardized Residuals GBM",
        names=c("M", "None", "F","Overlap"))
lines(c(0,46),c(2.5,2.5))
lines(c(0,46),c(-2.5,-2.5))
legend("topleft", legend=c("Depleted","Enriched"), fill=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), bty="n")

tnnGsx=table(NNdiaov[NNdiaov[,2]!=0,2]<0,sxn_nn[NNdiaov[,2]!=0])
chsqtnnGsxn=chisq.posthoc.test(tnnGsx, simulate.p.value=T)
chsqtnnGsxn=as.matrix(apply(chsqtnnGsxn[c(1,3),3:6],2,as.numeric))
barplot(chsqtnnGsxn, beside=T ,col=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), 
        xlab="Gender Related Compartments", ylab="Standardized Residuals Meta",
        names=c("M", "None", "F","Overlap"))
lines(c(0,46),c(2.5,2.5))
lines(c(0,46),c(-2.5,-2.5))
legend("topleft", legend=c("Depleted","Enriched"), fill=c(rgb(0.8,0.8,0.8,0.5),rgb(0,0,0,1)), bty="n")

posclqsx=lapply(posclq, function(l){
  x=sxn[names(l)]
  cut(x, c(-100,-0.5,0.5,100), labels=F)
})          

xM=t(sapply(posclqsx, function(x) c(length(x),sum(x==1)/length(x))))
aggregate(xM[,2], by=list(xM[,1]), "mean")
xF=t(sapply(posclqsx, function(x) c(length(x),sum(x==3)/length(x))))
aggregate(xF[,2], by=list(xF[,1]), "mean")

ctsx=t(sapply(ctspep,function(l) {
  tx=rep(0,3)
  names(tx)=c("-1","0","1")
  t0=table(sxn_n[l])
  tx[names(t0)]=t0
  return(tx)
}))
tclsx=table(calls=as.factor(pep %in% clcalls), sxn_n)
chisq.posthoc.test(as.table(ctsx), simulate.p.value=T)


# Cor. Sex - NN by clique -------------------------------------------------

mx=vb10-min(vb10)+1

ij=cut(seq_along(clqGvs),40, labels = F)
cl=makeCluster(4)
clusterExport(cl,c("clqGvs","mx","sex","ij"), envir = environment())
sexclq=pbsapply(1:40, function(i){
  cqi=clqGvs[ij==i]
  xx=as.data.frame(t(sapply(cqi, function(cq){
    cq=names(cq)
    y=c(mx[cq,])
    x=as.factor(rep(sex,each=length(cq)))
    d=data.frame(Intensity=y, Sex=x)
    tt=wilcox.test(data=d, Intensity~Sex)
    z=aggregate(Intensity~Sex,data=d,median)
    c(Stat=log10(z$Intensity[1]/z$Intensity[2]), P=tt$p.value)
  })))
  return(xx)
},cl=cl)
stopCluster(cl)

sexclq=data.frame(Stat=unlist(sexclq["Stat",]),P=unlist(sexclq["P",]))

sexclq$P=p.adjust(sexclq$P, method = "BH")

posexclq=clqGvs[sexclq$P<0.05]
lclqposex=lengths(posexclq)
lclqsex=lengths(clqGvs)
tlclqposex=table(lclqposex)
tlclqsex=table(lclqsex)
ntclqposex=as.numeric(names(tlclqposex))
ntclqsex=as.numeric(names(tlclqsex))
tlclqposex=c(tlclqposex, rep(0,length(ntclqsex)-length(ntclqposex)))
names(tlclqposex)=ntclqsex
plot(tlclqposex/tlclqsex, ylim=range(tlclqposex/tlclqsex), xlab="", ylab="", las=2)
par(new=T)
plot(c(3,38),c(0.4,0.4), ty="l", ylim=range(tlclqposex/tlclqsex), xlab="Clique Size", ylab="Proportion of Cliques with Sex Dependent Reactivtiy", las=2)

x=chisq.posthoc.test(cbind(tlclqsex-tlclqposex,tlclqposex), simulate.p.value = T)
j=(1:(nrow(x) %/% 2))*2
i=which(x[j,3]<0.05)

pepsexposclq=unlist(lapply(posexclq, names))
ttpepcqsx=table(pepsexposclq)
pepclq=unlist(lapply(clqGvs, names))
ttpepcq=table(pepclq)


ppsxoclq=names(ttpepcqsx)


pepbyclqsx=sapply(ppsxoclq, function(p){
  x1=ttpepcq[p]
  x2=ttpepcqsx[p]
  x=cbind(All=as.numeric(x1), Sex=as.numeric(x2))
  return(x)
})
nm=colnames(pepbyclqsx)
pepbyclqsx=t(data.matrix(pepbyclqsx))
rownames(pepbyclqsx)=nm
pepbyclqsxrt=pepbyclqsx[,2]/pepbyclqsx[,1]

clqsxFM=sexclq$Stat[sexclq$P<0.05]
pepsexFposclq=unlist(lapply(posexclq[clqsxFM>0], names))
ttpepcqFsx=table(pepsexFposclq)
pepsexMposclq=unlist(lapply(posexclq[clqsxFM<0], names))
ttpepcqMsx=table(pepsexMposclq)
xFM=union(names(ttpepcqFsx),names(ttpepcqMsx))
x=cbind(xFM %in% names(ttpepcqFsx),xFM %in% names(ttpepcqMsx))
rownames(x)=xFM
xFM=x
x=rep(0,length(vGvs))
names(x)=vGvs
x[rownames(xFM)[xFM[,1]==T]]=1
x[rownames(xFM)[xFM[,2]==T]]=-1
x[rownames(xFM)[xFM[,1]&xFM[,2]]]=0
table(x)

Gvs=set_vertex_attr(Gvs, "Sexclqpep", value =x)
write.graph(Gvs, format = "graphml", file="Gvs.graphml")

vgvsSexclq=x
vsxlq=as.factor(x)
boxplot(log10(dGvs)~vsxlq, notch=T, names=c("M","none","F"), ylab="Log Degree", xlab="Sex Differential Expression")
aovDsxclq=aov(log10(dGvs+0.5)~vsxlq)
summary(aovDsxclq)
TukeyHSD(aovDsxclq)

mnvs=rowMeans(vs10)
boxplot(mnvs~vsxlq, notch=T, names=c("M","none","F"), ylab="Mean", xlab="Sex Differential Expression")
aovMsxclq=aov(mnvs~vsxlq)
summary(aovMsxclq)
TukeyHSD(aovMsxclq)

ij=cut(rank(mnvs), 5, labels=F)
for (i in 1:5) {
  print(c(i,range(mnvs[ij==i])))
}
rtclqsx=as.table(aggregate(vsxlq, by=list(ij), "table")[[2]])
rownames(rtclqsx)=1:5
chsqclsx=chisq.posthoc.test(rtclqsx,simulate.p.value=T)
barplot(as.matrix(chsqclsx[chsqclsx[,2]=="Residuals",3:5]), ylab="Z Score", 
        xlab="Sex Differential Expression", beside=T, names=c("F","None","M"))
lines(c(0,18),c(3,3))
lines(c(0,18),c(-3,-3))
legend("topright",legend=paste("Mean Quintile ", 1:5, sep=""), 
       fill=c(rgb(0,0,0,0.875),rgb(0,0,0,0.75),rgb(0,0,0,0.5),rgb(0,0,0,0.25),rgb(0,0,0,0.125)),
       bty="n")


x=chsqclsx[chsqclsx[,2]=="p values",3:5]
x=apply(x,1,function(l) any(l<0.05))
chsqcls=chsqclsx[chsqclsx[,2]=="Residuals",][x,]


# Cor. AB0 -NN by clique --------------------------------------------------



ij=cut(seq_along(clqGvs),40, labels = F)
cl=makeCluster(4)
clusterExport(cl,c("clqGvs","mx","BlGr","ij"), envir = environment())
AB0clq=pbsapply(1:40, function(i){
  cqi=clqGvs[ij==i]
  xx=as.data.frame(t(sapply(cqi, function(cq){
    cq=names(cq)
    y=c(mx[cq,])
    xA=as.factor(rep(BlGr$A*1,each=length(cq)))
    xB=as.factor(rep(BlGr$B*1,each=length(cq)))
    x0=as.factor(rep(BlGr$N*1,each=length(cq)))
    d=data.frame(Intensity=y, A=xA, B=xB, N=x0)
    ttA=wilcox.test(data=d, Intensity~A)
    ttB=wilcox.test(data=d, Intensity~B)
    tt0=wilcox.test(data=d, Intensity~N)
    zA=aggregate(Intensity~A,data=d,median)
    zB=aggregate(Intensity~B,data=d,median)
    z0=aggregate(Intensity~N,data=d,median)
    resA=c(Stat=log10(zA$Intensity[2]/zA$Intensity[1]), P=ttA$p.value)
    resB=c(Stat=log10(zB$Intensity[2]/zB$Intensity[1]), P=ttB$p.value)
    res0=c(Stat=log10(z0$Intensity[2]/z0$Intensity[1]), P=tt0$p.value)
    return(data.frame(A=resA,B=resB, N=res0))
  })))
  return(xx)
},cl=cl)
stopCluster(cl)


xA=unlist(AB0clq[1,], recursive = F)
xB=unlist(AB0clq[2,], recursive = F)
x0=unlist(AB0clq[3,], recursive = F)
xA=t(sapply(xA, function(l) data.frame(Eff=l[1], p=l[2])))
xB=t(sapply(xB, function(l) data.frame(Eff=l[1], p=l[2])))
x0=t(sapply(x0, function(l) data.frame(Eff=l[1], p=l[2])))
AB0=cbind(xA,xB,x0)
colnames(AB0)=c("Aeff","Ap","Beff","Bp","Neff","Np")
rm(xA,xB,x0)
jA=p.adjust(AB0[,2])<0.05
jB=p.adjust(AB0[,4])<0.05
j0=p.adjust(AB0[,6])<0.05
jAE=AB0[,1]>0
jBE=AB0[,3]>0
j0E=AB0[,5]>0
ppAp=unique(unlist(lapply(clqGvs[jA&jAE], names)))
ppAn=unique(unlist(lapply(clqGvs[jA&(!jAE)], names)))
ppBp=unique(unlist(lapply(clqGvs[jB&jBE], names)))
ppBn=unique(unlist(lapply(clqGvs[jB&(!jBE)], names)))
pp0p=unique(unlist(lapply(clqGvs[j0&j0E], names)))
pp0n=unique(unlist(lapply(clqGvs[j0&(!j0E)], names)))
xAp=pep %in% ppAp
xAn=pep %in% ppAn
xBn=pep %in% ppBn
x0n=pep %in% pp0n
xAll=rowSums(cbind(xAp,2*xAn,4*xBn,8*x0n))
Gvs=set_vertex_attr(Gvs, "AB0clqpep", value =xAll)
write.graph(Gvs, format = "graphml", file="Gvs.graphml")


# ABO NN mean--------------------------------------------------------------


AB0mNN=t(sapply(vGvs, function(p){
    pp=names(ego(Gvs,1, p)[[1]])
    x=vs10[pp,]
    print(which(vGvs==p))
    if (!is.null(nrow(x))) n=nrow(x) else n=1
    Y=c(x)
    # xF=apply(cbind(BlGr[,1:2]*1,sex),1,paste,collapse="_")
    # xF=as.factor(rep(xF,each=n))
    pat=as.factor(rep(1:21,each=n))
    xA=as.factor(rep(BlGr$A,each=n))
    xB=as.factor(rep(BlGr$B,each=n))
    xS=as.factor(rep(sex,each=n))
    xS=summary(aov(Y~xS/pat*xA/pat*xB/pat)) 
    xA=summary(aov(Y~xA/pat))
    xB=summary(aov(Y~xB/pat))
    

    if(n>1) {y=c(F_A=as.numeric(xA[[1]][1,4]),
                 F_B=as.numeric(xB[[1]][1,4]),
                 F_S=as.numeric(xS[[1]][1,4]),
                 PadjA=as.numeric(xA[[1]][1,5]),
                 PadjB=as.numeric(xB[[1]][1,5]),
                 PadjS=as.numeric(xS[[1]][1,5]))
              return(y)
    }
    else {
      return(rep(1, 6))
    }
   
}))

mx=vb10-min(vb10)+1
Ltb=data.frame(Val=as.numeric(c(mx)), Sex=as.factor(rep(sex, each=nrow(vs10))), 
               A=as.factor(rep(BlGr$A*1+1,each=nrow(vs10))), 
               B=as.factor(rep(BlGr$B*1+1,each=nrow(vs10))))

x=aov(data=Ltb, Val~Sex+A*B)

AB0mNN[,4:6]=apply(AB0mNN[,4:6],2,p.adjust,method="BH",n=3*length(AB0mNN[,4]))
cor(log10(AB0mNN[AB0mNN[,4]<0.05&AB0mNN[,6]<0.05,c(1,3)]))
cor(log10(AB0mNN[AB0mNN[,5]<0.05&AB0mNN[,6]<0.05,c(2,3)]))
cor(as.numeric(as.factor(BlGr$A)),as.numeric(as.factor(sex)), method="spearman")
cor(as.numeric(as.factor(BlGr$B)),as.numeric(as.factor(sex)), method="spearman")
chisq.test(table(BlGr$A,sex), simulate.p.value = T)

chisq.test(table(AB0mNN[,4]<0.05,AB0mNN[,6]<0.05), simulate.p.value = T)

AB0sNN=t(sapply(vGvs, function(p){
  x=vs10[names(ego(Gvs,1, p)[[1]]),]
  print(which(vGvs==p))
  if (!is.null(nrow(x))) n=nrow(x) else n=1
  Y=c(x)
  xA=rep(BlGr$A,each=n)
  xB=rep(BlGr$B,each=n)
  xS=rep(sex,each=n)
  xA=mean(Y[xA])-mean(Y[!xA])
  xB=mean(Y[xB])-mean(Y[!xB])
  xS=mean(Y[xS=="F"])-mean(Y[xS=="M"])
  return(c(sign(xA),sign(xB), sign(xS)))
}))

AB0mNNg=AB0mNN[,1:3]
for (i in 1:3) {
  AB0mNNg[,i]=AB0mNNg[,i]*(AB0mNN[,i+3]<0.001)
}
for (i in 1:3) {
  AB0mNNg[,i]=AB0mNNg[,i]*AB0sNN[,i]
}

xa=cut((AB0mNN[,4]<0.001)*AB0sNN[,1], c(-2,-1e-4,1e-4,2), labels = c("Neg","None","Pos"))
xb=cut((AB0mNN[,5]<0.001)*AB0sNN[,2], c(-2,-1e-4,1e-4,23), labels = c("Neg","None","Pos"))
table(A=xa,B=xb)

Gvs=set_vertex_attr(Gvs, "AB0_A", value =AB0mNNg[,1])
Gvs=set_vertex_attr(Gvs, "AB0_B", value =AB0mNNg[,2])
Gvs=set_vertex_attr(Gvs, "AB0_S", value =AB0mNNg[,3])

write.graph(Gvs, format = "graphml", file="Gvs.graphml")

chisq.test(table(sex,BlGr$A), simulate.p.value = T)
chisq.test(table(sex,BlGr$B), simulate.p.value = T)
chisq.test(table(sex,BlGr$N), simulate.p.value = T)
table(sex,BlGr$N)
table(sex,BlGr$A)
table(sex,BlGr$B)

AB0SxKw=lapply(1:14, function(cl){
  j=names(GvcL)[GvcL==cl]
  n=length(j)
  Y=c(vs10[j,])
  xF=apply(cbind(BlGr*1,sex),1,paste,collapse="_")
  xF=as.factor(rep(xF,each=n))
  #boxplot(Y~xF, notch=T, main=cl, las=2, xlab="")
  x=dunn.test(Y,xF,method="bh")
  y=cbind(x$Z,x$P.adjusted)
  rownames(y)=x$comparisons
  colnames(y)=c("Effect","P.adjusted")
  y=y[y[,2]<0.05,]
  return(list(y[order(y[,1]),]))
})

# correlation of diagnosis associated reactivities

nnG=cut(NNdiaov[,1], c(-2,-1e-2,1e-2,2),labels=c("Gneg","Gnone","Gpos"))
nnM=cut(NNdiaov[,2], c(-2,-1e-2,1e-2,2),labels=c("Mneg","Mnone","Mpos"))
nnA=cut(AB0mNNg[,1], c(-1e3,-1e-1,1e-1,1e3),labels=c("Aneg","Anone","Apos"))
nnB=cut(AB0mNNg[,2], c(-1e3,-1e-1,1e-1,1e3),labels=c("Bneg","Bnone","Bpos"))

sink(file="diagXAB0.txt")
table(nnG,nnA)
chisq.test(table(nnG,nnA))
chisq.posthoc.test(table(nnG,nnA))
table(nnG,nnB)
chisq.test(table(nnG,nnB))
chisq.posthoc.test(table(nnG,nnB))
table(nnM,nnA)
chisq.test(table(nnM,nnA))
chisq.posthoc.test(table(nnM,nnA))
table(nnM,nnB)
chisq.test(table(nnM,nnB))
chisq.posthoc.test(table(nnM,nnB))
sink()

chisq.posthoc.test(table(x,pep %in% callsclqf), simulate.p.value=T)
chisq.posthoc.test(table(y,pep %in% callsclqf), simulate.p.value=T)
chisq.posthoc.test(table(AB0mNNg[,1]!=0,pep %in% callsclqf), simulate.p.value=T)
chisq.posthoc.test(table(AB0mNNg[,2]!=0,pep %in% callsclqf), simulate.p.value=T) # *
chisq.posthoc.test(table(AB0mNNg[,1]>0,pep %in% callsclqf), simulate.p.value=T)
chisq.posthoc.test(table(AB0mNNg[,1]<0,pep %in% callsclqf), simulate.p.value=T)
chisq.posthoc.test(table(AB0mNNg[,2]>0,pep %in% callsclqf), simulate.p.value=T)
chisq.posthoc.test(table(AB0mNNg[,2]<0,pep %in% callsclqf), simulate.p.value=T)
chisq.posthoc.test(table(apply(AB0mNNg[,1:2],1,function(l) any(l!=0)),pep %in% callsclqf), simulate.p.value=T)

mx=vb10-min(vb10)+1
AB0loNN=t(sapply(vGvs, function(p){
  pp=names(ego(Gvs,1, p)[[1]])
  x=mx[pp,]
  print(which(vGvs==p))
  if (!is.null(nrow(x))) {
    n=nrow(x)
    Y=colMeans(x)
  }
  else {
    n=1
    Y=x
  }
  pat=as.factor(1:21)
  xA=as.factor(BlGr$A)
  xB=as.factor(BlGr$B)
  xS=as.factor(sex)
  xA=diff(log2(aggregate(Y, by=list(xA),"mean")[,2]))
  xB=diff(log2(aggregate(Y, by=list(xB),"mean")[,2]))
  xS=diff(log2(aggregate(Y, by=list(xS),"mean")[,2]))
  return(c(A=xA, B=xB, S=xS))
}))

plot(AB0loNN[,3]~AB0loNN[,1])
plot(lm(AB0loNN[,3]~AB0loNN[,1])$residuals~AB0loNN[,2])

nnq10=table(cut(AB0loNN[,1],c(-10,-0.152,0.152,10),labels = c("Aneg","Anone","Apos")),
            cut(AB0loNN[,2],c(-10,-0.152,0.152,10),labels = c("Bneg","Bnone","Bpos")))
chisq.posthoc.test(nnq10, simulate.p.value=T)

SnoA=lm(AB0loNN[,3]~AB0loNN[,1]*AB0loNN[,2])
summary(SnoA)
SnoAr=SnoA$residuals
SnoArefd=ecdf(SnoAr)

si=quantile(SnoAr,c(0.1,0.9)) #c(0.01,0.05,0.1,0.5,0.9,0.95,0.99)
par(par0)
par(mai=c(1.5,1,1,0.5), mgp=c(4.2,1,0))
boxplot(rowMeans(vs10)~cut(SnoAr,c(-0.3,si,0.3)), 
       notch=T, xlab="Quantiles of the sex ratios (F/M)", ylab="Mean intensity", las=2 )#labels=c("<0.01","0.01-0.05","0.05-0.1","0.1-0.5","0.5-0.9","0.9-0.95","0.95-0.99",">0.99"))
par(par0)
x=aov(rowMeans(vb10)~cut(SnoAr,c(-0.3,si,0.3)))
TukeyHSD(x)
table(cut(SnoAr,c(-0.3,si[c(3,5)],0.3)))

Sdist=t(sapply(vGvs, function(p){
  pp=names(ego(Gvs,1, p)[[1]])
  x=mx[pp,]
  print(which(vGvs==p))
  if (!is.null(nrow(x))) {
    n=nrow(x)
    Y=colMeans(x)
  }
  else {
    n=1
    Y=x
  }
  xA=as.factor(BlGr$A)
  xB=as.factor(BlGr$B)
  xS=as.factor(sex)
  c(aggregate(Y, by=list(xS), "mean")[,2],
        aggregate(Y, by=list(xA), "mean")[,2],
        aggregate(Y, by=list(xB), "mean")[,2])
}))
colnames(Sdist)=c("Fm","M","negA","posA","negB","posB")
Sdist=as.data.frame(Sdist)
x=log2(Sdist$posA/Sdist$negA)
y=log2(Sdist$posB/Sdist$negB)
Sdist=data.frame(Sdist, Alo=x, Blo=y)
SdistF=lm(data=Sdist, Fm~posA+negA+posB+negB)$residuals
SdistM=lm(data=Sdist, M~posA+negA+posB+negB)$residuals
hF=hist(SdistF, xlim=range(c(SdistF,SdistM)), col=rgb(1,0,0,0.5), plot = F)
par(new=T)
hM=hist(SdistM, xlim=range(c(SdistF,SdistM)), col=rgb(0,0,1,0.5), plot = F)
plot(SdistF,SdistM)
plot(hF$mids, log10(hF$counts), ty="h", lwd=7, range(c(hF$mids,hM$mids)), 
     col=4, ylim=c(0,3.5), xlab="", ylab="")
par(new=T)
plot(hM$mids+0.015, log10(hM$counts), ty="h", lwd=7, range(c(hF$mids,hM$mids)), 
     col=2, ylim=c(0,3.5),
     xlab="Reactivity intensity", ylab="Number of neighborhoods")
legend("topright",legend = c("Women","Men"),fill=c(2,4),bty="n", cex=0.75)
par(new=F)

summary(lm(NNdiaov[,1]~SnoAr))
summary(lm(NNdiaov[,2]~SnoAr))
summary(lm(NNdiaov[,3]~SnoAr))

# Tumor associated reactivities were expressed stronger in women
# This effect was stronger for metastases.

summary(lm(NNdiaov[,1]~SnoAr*NNdiaov[,2]))

Gvs=set_vertex_attr(Gvs, "AB0loS", value =SnoAr)
Gvs=set_vertex_attr(Gvs, "AB0loA", value =AB0loNN[,1])
Gvs=set_vertex_attr(Gvs, "AB0loB", value =AB0loNN[,2])
write.graph(Gvs, format = "graphml", file="Gvs.graphml")

mnABS=sapply(1:10000, function(i){
  p=sample(pp,1)
  j=sample(21,7)
  a=aggregate(mx[p,], by=list((1:21) %in% j), "mean")
  diff(log2(a$x))
})
rA=range(mnABS)*0.5
mnBBS=sapply(1:10000, function(i){
  p=sample(pp,1)
  j=sample(21,4)
  a=aggregate(mx[p,], by=list((1:21) %in% j), "mean")
  diff(log2(a$x))
})
rB=range(mnBBS)*0.5

coefA=c(45,45,-45/(rA))
clrsA=1*(AB0loNN[,1]>0&-log10(AB0mNN[,4])>(coefA[2]+coefA[4]*AB0loNN[,1])&AB0mNN[,4]<0.05)+
     3*(AB0loNN[,1]<0&-log10(AB0mNN[,4])>(coefA[1]+coefA[3]*AB0loNN[,1])&AB0mNN[,4]<0.05)+1
plot(AB0loNN[,1],-log10(AB0mNN[,4]), ylim=c(0,70), pch=16, cex=0.35, col=clrsA, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="A") 
plot(AB0loNN[,1],-log10(AB0mNN[,4]), ylim=c(0,70), pch=16, cex=0.35, col=(abs(AB0loNN[,1])>0.152&AB0mNNg[,1]!=0)*1+1, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="A") 
plot(AB0loNN[,1],-log10(AB0mNN[,4]), ylim=c(0,70), pch=16, cex=0.2, col=(AB0mNN[,4]<0.001)*(1+2*(AB0mNNg[,1]<0))+1, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="Blood Group Antigen A") 

coefB=c(45,45,-45/(rB))
clrsB=1*(AB0loNN[,2]>0&-log10(AB0mNN[,5])>(coefB[2]+coefB[4]*AB0loNN[,2])&AB0mNN[,5]<0.05)+
     3*(AB0loNN[,2]<0&-log10(AB0mNN[,5])>(coefB[1]+coefB[3]*AB0loNN[,2])&AB0mNN[,5]<0.05)+1
plot(AB0loNN[,2],-log10(AB0mNN[,5]), ylim=c(0,70), pch=16, cex=0.35, col=clrsB, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="B")
plot(AB0loNN[,2],-log10(AB0mNN[,5]), ylim=c(0,70), pch=16, cex=0.35, col=(abs(AB0loNN[,2])>0.152&AB0mNNg[,2]!=0)*1+1, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="A") 
plot(AB0loNN[,2],-log10(AB0mNN[,5]), ylim=c(0,70), pch=16, cex=0.2, col=(AB0mNN[,5]<0.001)*(1+2*(AB0mNNg[,2]<0))+1, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="Blood Group Antigen B") 

clrsNNG=(NNdiaov[,1]<0)*2+(NNdiaov[,1]>0)*4
plot(AB0loNN[,1],-log10(AB0mNN[,4]), ylim=c(0,70), pch=16, cex=0.35, col=clrsNNG, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="A") 
plot(AB0loNN[,2],-log10(AB0mNN[,5]), ylim=c(0,70), pch=16, cex=0.35, col=clrsNNG, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="B")
clrsNNM=(NNdiaov[,2]<0)*2+(NNdiaov[,2]>0)*4
plot(AB0loNN[,1],-log10(AB0mNN[,4]), ylim=c(0,70), pch=16, cex=0.35, col=clrsNNM, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="A") 
plot(AB0loNN[,2],-log10(AB0mNN[,5]), ylim=c(0,70), pch=16, cex=0.35, col=clrsNNM, xlab="Log2 Fold Change", ylab="-Log10(p)", main ="B")


chisq.posthoc.test(table(cut(NNdiaov[,1], c(-100,-5e-1,5e-1,100), labels=F),cut(AB0loNN[,1], c(-100,-1.5e-1,1.5e-1,100), labels=F)), simulate.p.value=T)
chisq.posthoc.test(table(cut(NNdiaov[,2], c(-100,-5e-1,5e-1,100), labels=F),cut(AB0loNN[,1], c(-100,-1.5e-1,1.5e-1,100), labels=F)), simulate.p.value=T)
chisq.posthoc.test(table(cut(NNdiaov[,1], c(-100,-5e-1,5e-1,100), labels=F),cut(AB0loNN[,2], c(-100,-1.5e-1,1.5e-1,100), labels=F)), simulate.p.value=T)
chisq.posthoc.test(table(cut(NNdiaov[,2], c(-100,-5e-1,5e-1,100), labels=F),cut(AB0loNN[,2], c(-100,-1.5e-1,1.5e-1,100), labels=F)), simulate.p.value=T)

plot(table(cut(NNdiaov[,1], breaks=c(-100,-5e-1,5e-1,100), labels=c("Lost", "None", "Gained")),
           cut(AB0loNN[,1], breaks=c(-100,-1.5e-1,1.5e-1,100), labels=c("Lost", "None", "Gained"))), 
     main="Gliobalstoma by Blood Group A", xlab="Blood Group A", ylab="GBM") 

2^range(AB0loNN[AB0mNN[,1]<0.001&AB0loNN[,1]<0,1])
2^range(AB0loNN[AB0mNN[,1]<0.001&AB0loNN[,1]>0,1])
2^range(AB0loNN[AB0mNN[,2]<0.001&AB0loNN[,2]>0,2])

# Blood group cross-reactivtiy effects spread beyond cross-reactivity 
# boundaries - the farthest vertices of the induced subgraphs of 
# the A positive, B positive and A negative domains remain disconnected 
# up to depth 3 (2 for B negative)

ij=AB0mNNg[,1]>0
gvsA=induced.subgraph(Gvs, vGvs[ij])
diameter(gvsA, weights=rep(1,length(E(gvsA))))
x=names(farthest_vertices(gvsA)$vertices)
gpoles=names(unlist(ego(Gvs,3,x)))
x=induced.subgraph(Gvs, gpoles)
components(x)

ij=AB0mNNg[,2]>0
gvsA=induced.subgraph(Gvs, vGvs[ij])
diameter(gvsA, weights=rep(1,length(E(gvsA))))
x=names(farthest_vertices(gvsA)$vertices)
gpoles=names(unlist(ego(Gvs,2,x)))
x=induced.subgraph(Gvs, gpoles)
components(x)

# Modularities ------------------------------------------------------------

### sex

sxj=cut(sxnn, c(-1e16, -1e-16,1e-16, 1e16), labels = c(-1,0,1))
sexmod=Zmod(Gvs, sxj)
CV=rowSds(vb010[names(V(Gvs)),])/rowMeans(vb010)[names(V(Gvs))]
CV=cut(CV, quantile(CV, c(0,0.33,0.67,1))-c(0.001,0,0,-0.001), labels=F)
cvmod=Zmod(Gvs, CV)
mn=cut(rowMeans(vb10), quantile(rowMeans(vb10), c(0,0.33,0.67,1))-c(0.001,0,0,-0.001), labels=F)
mnmod=Zmod(Gvs,mn)

ij=cut(rank(mnvs), 5, labels=F)
sxmnmod=sapply(1:5,function(i){
  jj=vsxlq %in% c(-1,1)
  j=(ij==i)[jj]
  vf=factor(vsxlq,levels=c(-1,1))
  vf=vf[!is.na(vf)]
  gi=induced.subgraph(Gvs, vGvs[jj][j])
  Zmod(gi, vf[j])
})

sxmnass=sapply(1:5,function(i){
  jj=vsxlq %in% c(-1,1)
  j=(ij==i)[jj]
  vf=factor(vsxlq,levels=c(-1,1))
  vf=vf[!is.na(vf)]
  gi=induced.subgraph(Gvs, vGvs[jj][j])
  Zass(gi, vf[j])
})

sexclqmod=Zmod(Gvs, as.factor(vgvsSexclq))

### degree
degtr=cut(dGvs, quantile(dGvs, c(0,0.33,0.67,1))-c(0.001,0,0,-0.001), labels=F)
degass=Zass(Gvs, degtr)
degmod=Zmod(Gvs, degtr)

# AB0

Y=AB0mNNg!=0
AB0assA=Zass(Gvs, Y[,1])            
AB0assB=Zass(Gvs, Y[,2])
AB0modA=Zmod(Gvs, Y[,1]*1+1)            
AB0modB=Zmod(Gvs, Y[,2]*1+1)
