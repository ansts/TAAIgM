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
require(Rfast)
require(lsmeans)
require(nlme)
require(heatmap3)
require(vioplot)
require(ggplot2)

require(irr)
require(pROC)
require(uwot)
require(RColorBraewer)
require(qualV)
require(infotheo)

require(rgl)
require(corrplot)
require(Peptides)
require(Biostrings)
require(factoextra)
require(FactoMineR)
require(riverplot)
require(class)



load("IgJtrim")
viral=c(2,4,9:18,20,21,26:35,41,42,45,46,48:50,52,55,60)
vir=pepinfo$peptide[pepinfo$protein %in% pepi[viral]]
palette(c("black","red","green","blue","yellow","cyan","magenta","gray75","gray35","dodgerblue2","firebrick4","darkolivegreen","lightcyan2","hotpink4","ivory3","gray90"))


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

base=c(vb10)[which(rank(vb10)==2)]-min(vb10)
mx=vb10-min(vb10)+base

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


# Simulation of profile diversity -----------------------------------------

prfls=t(sapply(pep,function(p){
  x=names(ego(Gvs,1,p)[[1]])
  m=mx[x,]
  if (length(x)>1) (colMeans(m)>mean(m))*1 else (m>mean(m))*1
}))

cl=makeCluster(14)
clusterExport(cl,"Groupings")
ps=t(pbapply(prfls,1,function(l){
  apply(Groupings, 2,function(g){
    chisq.test(l,g, simulate.p.value=T)$p.value
  })
},cl=cl))
stopCluster(cl)

ps=apply(ps,2,p.adjust)

# Cliques -----------------------------------------------------------------

clqGvs=max_cliques(Gvs, min=3) 
tclql=table(lengths(clqGvs)) 
vgrep=Vectorize(grep, vectorize.args = "pattern")
ls=as.numeric(names(tclql))

clqord=clqGvs[order(lengths(clqGvs))]
x=list()
xp=c()
for (i in seq_along(clqord)) {
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

hist(log10(dGvs[unique(unlist(uniclq))]+0.5), xlim=range(log10(dGvs+0.5)))
par(new=T)
hist(log10(dGvs+0.5), col=2, xlim=range(log10(dGvs+0.5)))

bigclq=unlist(lapply(clqord[lengths(clqord)>20], function(clq) names(clq)))
tbigclq=table(bigclq)
smallclq=unlist(lapply(clqord[lengths(clqord)<5], function(clq) names(clq)))
tsmallclq=table(smallclq)
length(intersect(names(tbigclq), names(tsmallclq)))
length(tbigclq)
length(tsmallclq)
sum(names(tbigclq) %in% vir)
sum(names(tsmallclq) %in% vir)
tbigclq[names(tbigclq) %in% vir]
chisq.posthoc.test(table(names(tbigclq) %in% vir,tbigclq>100), simulate.p.value=T)
X=lapply(clqord, names)
X=melt(X)
X=cbind(X, idio=q7ppp[X$value])
i7=aggregate(X$idio, by=list(X$L1), "mean")

bsi7=sapply(1:38, function(l){
  summary(sapply(1:5000, function(i){
    p=sample(pep,l)
    mean(q7ppp[p])
  }))
})
y=log10(i7$x+0.5)
vioplot(y~lengths(clqord)[i7$Group.1], ylim=range(y), xlab = "Clique size", 
        ylab="Log10(Mean # idiotope hits + 0.5) ")
par(new=T)
plot(3:38,log10(bsi7[2,3:38]+0.5), xlim=c(2.5,38.5), ylim=range(y), ty="l", 
     lwd=2, col=3, xlab="", ylab="", xaxt="n")
par(new=T)
plot(3:38,log10(bsi7[3,3:38]+0.5), xlim=c(2.5,38.5), ylim=range(y), ty="l", 
     lwd=2, col=2, xlab="", ylab="", xaxt="n")
par(new=T)
plot(3:38,log10(bsi7[5,3:38]+0.5), xlim=c(2.5,38.5), ylim=range(y), ty="l", 
     lwd=2, col=3, xlab="", ylab="", xaxt="n")

log10(median(q7ppp[names(tbigclq)])+0.5)
log10(median(q7ppp[names(tsmallclq)])+0.5)
chisq.posthoc.test(table(as.numeric(tbigclq)>20,q7ppp[names(tbigclq)]>0))
chisq.posthoc.test(table(as.numeric(tsmallclq)>3,q7ppp[names(tsmallclq)]>0))

pepbycl=table(unlist(lapply(clqord, names)))
pepbyclpo=table(unlist(lapply(posclqall, names)))

cl=makeCluster(18)
clusterExport(cl,c("vs10", "dgnf", "clucri","connfix","dunnfix","bhgamfix","BHgamma"))
j=sample(seq_along(clqord))
clqvirclucr=pbsapply(clqord[j], function(cq){
  cq=names(cq)
  l=length(cq)
  if (l==3) {
    x=clucri(vs10[cq,], dgnf)
  }
  else {
    cm=combn(l,3)
    n=ncol(cm)
    if (n>100) j=sample(1:n,100) else j=1:n
    x=mean(apply(cm[,j],2,function(i){
      clucri(vs10[cq[i],], dgnf)
    }))
  }
},cl=cl)
stopCluster(cl)

X=lapply(clqord, names)
X=melt(X)
ij=table(X$L1)
X=cbind(X,L=ij[X$L1])
X=X[,-3]
X=aggregate(X$L.Freq, by=list(X$value),"max")
Xclucr=sapply(3:20, function(n){
  clucri(vs10[X$Group.1[X$x==n],], dgnf)
})
plot(Xclucr, ty="b")
y=cut(X$x,c(0,3,4,5,6,7,8,9,10,12,15,20,24,38), labels=F)

xvir=sapply(clqord[j], function(cq) {
  cq=names(cq)
  sum(cq %in% vir)/length(cq)
})

xvirn=sapply(clqord[j], function(cq) {
  cq=names(cq)
  sum(cq %in% vir)
})

clqvirclucr=cbind(Crit=clqvirclucr, L=lengths(clqord[j]), Vir=xvir, Virn=xvirn, Idiotype=)
boxplot(data=clqvirclucr, Crit~L, notch=T, cex=0.3)
xv=cut(log10(clqvirclucr[,3]+0.001),9, labels=FALSE)
x=cut(clqvirclucr[,1],10,labels=F)
vioplot(clqvirclucr[,3]~x)
plot3d(clqvirclucr[,1:3])
rglwidget()
plot3d(clqvirclucr[,c(1,2,4)])
rglwidget()

jj=cut(clqvirclucr[,3], 3, labels=F)
for (i in 1:3) {
  plot(clqvirclucr[jj==i,c(2,1)], pch=16, xlim=c(3,38), ylim=c(-10,6),col=rgb(0.5,0.5,0.5,0.5), main=i)
}

critqnt=quantile(clqvirclucr[,1], c(0.33,0.67))
x=aggregate(clqvirclucr[,3], by=list(clqvirclucr[,2]), "mean")
x33=aggregate(clqvirclucr[clqvirclucr[,1]<critqnt[1],3], by=list(clqvirclucr[clqvirclucr[,1]<critqnt[1],2]), "mean")
x67=aggregate(clqvirclucr[clqvirclucr[,1]>critqnt[2],3], by=list(clqvirclucr[clqvirclucr[,1]>critqnt[2],2]), "mean")

plot(x, xlab="Crossreactivity clique size", ylab="Mean proportion of viral epitopes", xlim=c(3,38), ylim=range(x[,2]), ty="b")
par(new=T)
plot(x33, xlab="", ylab="", xlim=c(3,38), ylim=range(x[,2]), col=4, ty="l")
par(new=T)
plot(x67, xlab="", ylab="", xlim=c(3,38), ylim=range(x[,2]), col=2, ty="l")
lines(c(3,38),c(mean(x[,2]),mean(x[,2])))
legend("bottomright", legend=c("All max cliques >2","1st Tertile clustering coeff.","3rd Tertile clustering coeff."), lwd=1, col=c(1,4,2), bty="n")

lmvirclq=lm(Crit~Vir*L, data = as.data.frame(clqvirclucr))
summary(lmvirclq)

x=aggregate(clqvirclucr[,3], by=list(clqvirclucr[,2]), "mean")
y=aggregate(clqvirclucr[,1], by=list(clqvirclucr[,2]), "mean")
xy=cbind(x,y$x)
colnames(xy)=c("L","V","C")
summary(lm(data=xy, C~V*L))
summary(lm(data=xy,V~L))

Gvsnovir=induced.subgraph(Gvs, vGvs[!(vGvs %in% vir)])
clqsnovir=max_cliques(Gvsnovir, min=3)

cl=makeCluster(15)
clusterExport(cl, c("Gvs","vGvs","vir"))
clusterEvalQ(cl, require(igraph))
clqGnovirBS=t(pbsapply(1:1000, function(i){
  print(i)
  tb0=rep(0,36)
  names(tb0)=3:38
  x=sample(vGvs,length(vir))
  g=induced_subgraph(Gvs, vGvs[!(vGvs %in% x)])
  clq=max_cliques(g, min=3)
  tb=table(lengths(clq))
  tb0[names(tb)]=tb
  return(tb0)
}, cl=cl))
stopCluster(cl)
clqGnovirBS=sweep(clqGnovirBS,1,rowsums(clqGnovirBS),"/")

write_graph(Gvsnovir, format = "graphml", file="Gvsnovir.graphml")
clqsnovirL=lengths(clqsnovir)

tbvir=table(clqsnovirL)/length(clqsnovirL)
boxplot(clqGnovirBS, ylim=range(tbvir), range=1, outline=F, xlab="Clique size", ylab="Probability")
par(new=T)
plot(tclql/sum(tclql), ty="l", col=3, ylim=range(tbvir), xlim=c(2.5,38.5), xlab="", ylab="")
par(new=T)
plot(tbvir, ty="l", col=2, xlim=c(2.5,38.5), ylim=range(tbvir), xlab="", ylab="")
legend("topright", legend=c("Entire graph","The graph without viral epitopes","Bootstrap of subgraph"), bty="n", lwd=c(1,1,10),  col=c(3,2,8), cex = 1)

virb=vir[vir %in% names(tbigclq)]
virs=vir[vir %in% names(tsmallclq)]
bigsmallclq=data.frame(`Log(N)`=c(log10(tbigclq[!(names(tbigclq) %in% vir)]), log10(tbigclq[virb]),
                    log10(tsmallclq[!(names(tsmallclq) %in% vir)]), log10(tsmallclq[virs])),
                  Size=c(rep("Big", length(tbigclq[!(names(tbigclq) %in% virb)])+ length(tbigclq[virb])),
                         rep("Small", length(tsmallclq[!(names(tsmallclq) %in% vir)])+length(tsmallclq[virs]))),
                  Self=c(rep("Self", length(tbigclq[!(names(tbigclq) %in% virb)])),
                         rep("Vir", length(tbigclq[virb])),
                         rep("Self",length(tsmallclq[!(names(tsmallclq) %in% vir)])),
                         rep("Vir",length(tsmallclq[virs]))))
boxplot(data=bigsmallclq, Log.N.~Self*Size, notch=T)
summary(lm(data=bigsmallclq, Log.N.~Self*Size))
pepsmallbig=data.frame(Pep=c(names(tbigclq[!(names(tbigclq) %in% vir)]),names(tbigclq[virb]),
                             names(tsmallclq[!(names(tsmallclq) %in% vir)]),names(tsmallclq[virs])),
                       Size=c(rep("Big", length(tbigclq[!(names(tbigclq) %in% virb)])+ length(tbigclq[virb])),
                              rep("Small", length(tsmallclq[!(names(tsmallclq) %in% vir)])+length(tsmallclq[virs]))),
                       Self=c(rep("Self", length(tbigclq[!(names(tbigclq) %in% virb)])),
                              rep("Vir", length(tbigclq[virb])),
                              rep("Self",length(tsmallclq[!(names(tsmallclq) %in% vir)])),
                              rep("Vir",length(tsmallclq[virs]))))
chisq.posthoc.test(table(pepsmallbig$Size,pepsmallbig$Self), simulate.p.value=T)
bigclqpep=names(tbigclq)
smallclqpep=names(tsmallclq)

# Entropies

bcqpp=sapply(bigclqpep,function(p) unlist(strsplit(p, split="")))
scqpp=sapply(smallclqpep,function(p) unlist(strsplit(p, split="")))
vbcqpp=sapply(virb,function(p) unlist(strsplit(p, split="")))
vscqpp=sapply(virs,function(p) unlist(strsplit(p, split="")))

bcqppent=sapply(bcqpp, function(p) seqient(seqdef(t(p)), norm=F, base = 2))
scqppent=sapply(scqpp, function(p) seqient(seqdef(t(p)), norm=F, base = 2))
vbcqppent=sapply(vbcqpp, function(p) seqient(seqdef(t(p)), norm=F, base = 2))
vscqppent=sapply(vscqpp, function(p) seqient(seqdef(t(p)), norm=F, base = 2))

ents=list(BigCliques=bcqppent,SmallCliques=scqppent,Viral_BClq=vbcqppent,Viral_SClq=vscqppent)
boxplot(ents, notch = T)
ents=melt(ents)
colnames(ents)=c("value","Epitope_Set")
ggplot(data=ents,aes(x = Epitope_Set, y = value))+
  geom_violin(trim=T,fill='#A0A090',color="darkgreen") + 
  theme_classic()+
  stat_summary(fun = "mean",
               geom = "crossbar",
               width=0.5,
               color = "black")+
  ylab("Entropy [bits]")
ents$Epitope_Set=as.factor(ents$Epitope_Set)
entslm=lme(data=ents, value~Epitope_Set, random = ~1|value)
resentslm=lsmeans(entslm, pairwise~Epitope_Set)
resentslm=summary(resentslm)

xs=sapply(smallclqpep,as.AAbin)
xb=sapply(bigclqpep,as.AAbin)
xsd=kdistance(xs, k=4)
xbd=kdistance(xb, k=4)

xsdc=cut(xsd, c(-1,0.05,0.1,0.2,0.4,0.5,0.6,0.7,1), labels=F)
xbdc=cut(xbd, c(-1,0.05,0.1,0.2,0.4,0.5,0.6,0.7,1), labels=F)
txdc=cbind(small=table(xsdc), big=table(xbdc))

chisq.posthoc.test(txdc, simulate.p.value = T)

plot(log10(sweep(txdc,2,colSums(txdc),"/")), xlim=c(-4,0), ylim=c(-4,0))



#################

x=names(pepiused[pepiused %in% pepi[16:17]])
xbcq=intersect(x,bigclqpep)
y=length(bigclqpep)
tbxclq=rbind(c(length(pep)-length(x),length(x)),c(y-length(xbcq),length(xbcq)))
chsqtxbclq=chisq.test(tbxclq, simulate.p.value=T)
chsqtxbclq$stdres
pepbyclqtb=data.frame(Prots=pepiused[pep], Nclq=pepbycl[pep], Vir=pep %in% vir)
pepbyclqtb=pepbyclqtb[,-2]
pepbyclqtb$Prots[is.na(pepbyclqtb$Prots)]="Linker"
pepbyclqtb$Nclq.Freq[is.na(pepbyclqtb$Nclq.Freq)]=0
x=cut(pepbyclqtb$Nclq.Freq, c(-1,5,21,1e6), labels=F)
pepbyclqtb=data.frame(pepbyclqtb,NclqLT=x)
rownames(pepbyclqtb)=pep
summary(lm(data=pepbyclqtb, Nclq.Freq~Vir))
chisq.posthoc.test(table(pepbyclqtb$Virs,pepbyclqtb$NclqLT.1), simulate.p.value=T)

y=aggregate(data=pepbyclqtb, pep~Prots*NclqLT, "length")
y=dcast(data=y, Prots~NclqLT)
y[is.na(y)]=0
rownames(y)=y$Prots
y=y[,-1]
yvir=y[pepi[viral],]
virus=pepi[viral]
virs=c("EBV","EBV","HPV","HPV","HPV","HPV","HPV","HPV","EBV","HTLV","HTLV","EBV","HBV","KSHV","HPV","HPV","HPV","HPV","KSHV","EBV","HPV","Hp","HTLV","HPV","HBV","HPV","HTLV","EBV","KSHV","HTLV")
names(virs)=virus
virtab=aggregate(yvir, by=list(virs[rownames(yvir)]),"sum")
rownames(virtab)=virtab$Group.1
virtab=virtab[,-1]
chisq.posthoc.test(virtab, simulate.p.value=T)
pepbyclqtb=data.frame(pepbyclqtb,Virs=virs[pepbyclqtb$Prots])
pepbyclqtb$Virs[is.na(pepbyclqtb$Virs)]="self"
pepbyclqtb$Virs=as.factor(pepbyclqtb$Virs)
pepbyclqtb$Virs=relevel(pepbyclqtb$Virs,"HTLV")
virlm=lm(data=pepbyclqtb, subset = Vir==T, Nclq.Freq~Virs)
summary(virlm)

seqtoFASTA(bigclqpep, "bigclq.txt")

eiGvs=eigen_centrality(Gvs, directed = F,weights=NA)
plot(log10(pepbyclqtb[vGvs,2]+0.5)+rnorm(length(dGvs),0,0.3), log10(eiGvs$vector+min(eiGvs$vector)/2), pch=16, cex=0.6,col=rgb(0,0,0,0.3))
x=log10(as.numeric(unlist(eiGvs$vector))+1e-11)
Gvs=set.vertex.attribute(Gvs, "EigenCentr", value=x)
write.graph(Gvs, format = "graphml", file="Gvs.graphml")

x=names(tbigclq)
y=names(tsmallclq)
clqsizeGrps=data.frame(melt(vs10),
                       apply(Groupings,2,function(clm) as.factor(rep(clm, each=nrow(vs10)))),
                       BigSmall=rep((pep %in% y)+(pep %in% x)*2,21), Viral=as.factor(1+(pep %in% vir)))
clqsizeGrps$BigSmall[clqsizeGrps$BigSmall>1]=5-clqsizeGrps$BigSmall[clqsizeGrps$BigSmall>1]
clqsizeGrps$BigSmall=as.factor(clqsizeGrps$BigSmall)
colnames(clqsizeGrps)[1:2]=c("Seq","Pat")

AclqszVir=lme(data=clqsizeGrps, value~A*Viral*BigSmall,random=~1|Pat)
resAClqsV=lsmeans(AclqszVir, pairwise~A*Viral|BigSmall)
resAClqsV=summary(resAClqsV)

X=aggregate(clqsizeGrps$value,by=list(clqsizeGrps[,2],clqsizeGrps[,4],clqsizeGrps[,5],
                                      clqsizeGrps[,6],clqsizeGrps[,7],clqsizeGrps[,8],
                                      clqsizeGrps[,9],clqsizeGrps[,10],clqsizeGrps[,11],
                                      clqsizeGrps[,12]), "mean")
colnames(X)[2:10]=colnames(clqsizeGrps)[4:12]
colnames(X)[11]="Reactivity"
for (i in 2:10) X[,i]=as.factor(X[,i])

#Fr0=as.factor(clqsizeGrps$BigSmall %in% c("2","3")+1)
Fr0=as.factor((clqsizeGrps$BigSmall==1)+1)
AclqszVir=lme(data=clqsizeGrps, value~A*Viral*Fr0,random=~1|Pat)
resAClqsV=lsmeans(AclqszVir, pairwise~A*Viral|Fr0)
resAClqsV=summary(resAClqsV)
BclqszVir=lme(data=clqsizeGrps, value~B*Viral*Fr0,random=~1|Pat)
resBClqsV=lsmeans(BclqszVir, pairwise~B*Viral|Fr0)
resBClqsV=summary(resBClqsV)

#Fr=as.factor(X$BigSmall %in% c("2","3")+1)
Fr=as.factor((X$BigSmall=="2")+1)
Groups=as.factor(apply(X[,c(2,10)],1,function(x){
  if (x[1]) A="A" else A="nonA"
  if (x[2]==2) V="Viral" else V="Self"
  paste(A,V,sep="_")
}))
Xi=cbind(X[Fr==2,],Groups=Groups[Fr==2])
ggplot(data=Xi,aes(x = Groups, y = Reactivity))+
  geom_violin(trim=T,fill='#A0A090',color="darkgreen") + 
  theme_classic()+
  stat_summary(fun = "mean",
               geom = "crossbar",
               width=0.5,
               color = "black")+
  ggtitle("Big cliques")+
  ylab("Reactivity (a.u.)")


# Big Cliques Viral Subgraph ----------------------------------------------

x=unique(unlist(sapply(ego(Gvs, 1, names(tbigclq)),names)))
Gbclvir=induced_subgraph(Gvs, x)
vGbclvir=names(V(Gbclvir))
Gbclvir=set_vertex_attr(Gbclvir,"Viral",value=(vGbclvir %in% vir)+1)
Gbclvir=set_vertex_attr(Gbclvir,"BigCl_core",value=(vGbclvir %in% names(tbigclq))+1)
write.graph(Gbclvir, format = "graphml", file="Gbclvir.graphml")


# Cliques 2 ---------------------------------------------------------------

clqGvs10=max_cliques(Gvs, min=3, max=10)
mx=vb10-min(vb10)+1

diaclq=as.data.frame(t(sapply(clqGvs10, function(cq){
  cq=names(cq)
  y=colMeans(mx[cq,])
  c(Stat=kruskal.test(y~dgnf)$statistic, P=kruskal.test(y~dgnf)$p.value)
})))

diaclq$P=p.adjust(diaclq$P, method = "BH")


barplot(tclql, xlab= "Clique Size", ylab="N", yaxt="n", xaxt="n")
axis(1, cex.axis=0.75, at=(0:35)*1.2+0.7, labels=3:38, las=2)
axis(2, cex.axis=0.75)
barplot(tposclqall, xlab= "Clique Size", ylab="N", yaxt="n", xaxt="n")
axis(1, cex.axis=0.75, at=(0:35)*1.2+0.7, labels=3:38, las=2)
axis(2, cex.axis=0.75)

N=seq_along(mx[,1])
cl=makeCluster(15)
clusterExport(cl, c("mx","dgnf","clucri","connfix","dunnfix","bhgamfix","BHgamma","N"), envir = environment())
bsx=pblapply(3:10, function(n){
  sapply(1:10000,function(i){
    bscq=sample(N,n)
    clucri(mx[bscq,],dgnf)
  })
},cl=cl)
stopCluster(cl)
names(bsx)=3:10

cl=makeCluster(15)
names(clqGvs10)=seq_along(clqGvs10)
nY=names(clqGvs10)
nY=sample(nY)
clusterExport(cl, c("mx","dgnf","clucri","connfix","dunnfix","bhgamfix","BHgamma","clqGvs10"), envir = environment())
clucri3posclq=t(pbsapply(nY, function(ni){
    clq=clqGvs10[[ni]]
    x=names(clq)
    cc=clucri(mx[x,], dgnf)
    return(c(length(x),cc))
}, cl=cl))
stopCluster(cl)

crclcr=lapply(3:10,function(i) ecdf(bsx[[as.character(i)]]))
  
j=apply(clucri3posclq,1,function(l) {
  i=l[1]-2
  print(i)
  F=crclcr[[i]]
  F(l[2])>0.99
})

selposclq=clucri3posclq[j,]
colnames(selposclq)=c("N","Crit")
posclq=rownames(selposclq)
posclqpep=unique(unlist(sapply(clqGvs10[posclq],names)))
  
  #boxplot(clucri3posclq[,1]~clucri3posclq[,2], notch=T)
  
  dgnGf=dgnG&(sex=="F")
  dgnGm=dgnG&(sex=="M")
  dgnfm=dgnf
  dgnfm[dgnGf]=4

selFposclqp=Loolocal0(posclqpep,mx)
plotMDS(vs10[selFposclqp,], col=dgnf, pch=16)

length(intersect(selFposclqp,names(clqFTb)))/length(selFposclqp)
length(intersect(selFposclqp,callsclqf))/length(selFposclqp)

psclq=sapply(clqGvs10[posclq],names)
dmxposclq=sapply(psclq[order(lengths(psclq))], function(cl1){
            sapply(psclq[order(lengths(psclq))], function(cl2){
              as.numeric(length(intersect(cl1,cl2)))
            })
})

Gpcl=graph_from_adjacency_matrix(dmxposclq, mode = "undirected", weighted = TRUE, diag = FALSE)
vGpcl=names(V(Gpcl))
dGpcl=degree(Gpcl)
clqGpcl=max_cliques(Gpcl)
plot(table(lengths(clqGpcl)))
sz=lengths(psclq[vGpcl])
Groupings=cbind(Groupings, sex=(sex=="F"))
AttrGpcl=t(sapply(psclq[order(lengths(psclq))], function(cl){
  m=c(mx[cl,])
  n=length(cl)
  at=apply(Groupings,2,function(a){
    a=as.factor(rep(a, each=n))
    x=aggregate(m, by=list(a), "mean")$x
    x[2]/x[1]
  })
}))
AttrGpcl=scale(AttrGpcl)

for (j in 1:7) {
  at=AttrGpcl[,j]
  Gpcl=set.vertex.attribute(Gpcl, name=colnames(Groupings)[j], value=at)
}
Gpcl=set.vertex.attribute(Gpcl, "Size", value=sz)
write.graph(Gpcl, format = "graphml", file="Gpcl.graphml")



# Louvain clustering on Gpcl

GpclqLouv=cluster_louvain(Gpcl) 
GpclqLouv=aggregate(vGpcl,by=list(GpclqLouv$memberships[2,]),"list")$x
names(GpclqLouv)=seq_along(GpclqLouv)
x=melt(GpclqLouv)
xn=x$value
x=x$L1
names(x)=xn

glvmap=lapply(GpclqLouv, function(cl){
  px=clqclqpp$Pep[as.numeric(clqclqpp$CLq1) %in% as.numeric(cl)]
  table(GvcL[px])
})

Gpcl=set.vertex.attribute(Gpcl, "Louv", value=x[vGpcl])
write.graph(Gpcl, format = "graphml", file="Gpcl.graphml")

clqclq=melt(GpclqLouv)
j=duplicated(clqclq[,1])
clqclqpp=apply(clqclq,1,function(x){
  n=as.numeric(x[1])
  y0=psclq[[as.character(n)]]
  l=length(y0)
  x1=unlist(rep(x[1],l))
  x2=unlist(rep(x[2],l))
  y=data.frame(Pep=y0,CLq1=x1, Clq2=x2)
  print(list(x,y0,y))
  return(y)
})
x=c()
for (l in clqclqpp){
  x=rbind(x,l)
}
clqclqpp=x
clqclqp=clqclqpp[,-2]
clqclqp=unique(clqclqp)
clqclqpdf=cbind(clqclqp, vs10[clqclqp$Pep,])
clqclqpdf=melt(clqclqpdf)

colnames(clqclqpdf)[3]="Pat"
P=as.character(clqclqpdf$Pat)
clqclqpdf=data.frame(clqclqpdf, A=as.factor(Groupings[P,1]),
                     B=as.factor(Groupings[P,2]),
                     O=as.factor(Groupings[P,3]),
                     C=as.factor(Groupings[P,4]),
                     G=as.factor(Groupings[P,5]),
                     M=as.factor(Groupings[P,6]),
                     S=as.factor(Groupings[P,7]))
clqclqpdf$Clq2=as.factor(clqclqpdf$Clq2)

V=as.numeric(clqclqpdf$value)
Clq=as.factor(clqclqpdf$Clq2)
P=as.factor(clqclqpdf$Pat)
rres=lapply(5:11, function(i){
  print(i)
  Par=as.factor(clqclqpdf[,i])
  Clqclqlm=lme(V~Par*Clq, random=~1|P)
  res=lsmeans(Clqclqlm, pairwise~Par|Clq)
  reres=summary(res$contrasts)
  j=reres$p.value<0.1
  rres=data.frame(Contrast=as.character(reres$contrast[j]),Clique=as.character(reres$Clq[j]),Estimate=as.numeric(reres$estimate[j]),p.value=as.numeric(reres$p.value[j]))
  rres[order(rres$Estimate),]
})
names(rres)=colnames(clqclqpdf)[5:11]
sink(file="Xnew.txt")
print(rres)
sink()
ij=rbind(rep(5:7, each=4),rep(8:11,3))
rres3=apply(ij, 2,function(i){
  print(i)
  P1=as.factor(clqclqpdf[,i[1]])
  P2=as.factor(clqclqpdf[,i[2]])
  PxP=lme(V~P1*P2*Clq,random=~1|P)
  resPP=lsmeans(PxP, pairwise~P1|P2|Clq)
  resPP=summary(resPP)
  resPP=resPP$contrasts[resPP$contrasts$p.value<0.05,]
  resPP=resPP[order(resPP$P2,resPP$estimate),]
  colnames(resPP)[1]=paste(colnames(clqclqpdf)[i[1]],"contrast",sep="_")
  colnames(resPP)[2]=colnames(clqclqpdf)[i[2]]
  return(resPP)
})
sink(file="Clq2Louvxparams_3way_new.txt")
print(rres3)
sink()



# Stats of clq2 induced subgraph ------------------------------------------

clq2=rep(0,length(vGvs))
names(clq2)=vGvs
clq2[clqclqpp$Pep]=clqclqpp$Clq2
Gvs=set.vertex.attribute(Gvs, name = "Clq2new", value=clq2)
write.graph(Gvs, format = "graphml", file="Gvs.graphml")

chsqgvlclq2=chisq.posthoc.test(table(GvcL[vGvs],clq2), simulate.p.value=T)
ij=(1:(nrow(chsqgvlclq2)/2))*2
ji=ij-1
ii=which(chsqgvlclq2[ij,]<0.05,arr.ind = T)
i=chsqgvlclq2$Dimension[ij][ii[,1]]
j=colnames(chsqgvlclq2)[ii[,2]]
chsqPHres=data.frame(Cluster=i, Clique2Cl=j,Eff=as.numeric(chsqgvlclq2[ji,][ii]))
chsqPHres=chsqPHres[order(as.numeric(chsqPHres$Cluster), chsqPHres$Eff),]

Gvsclq2=induced.subgraph(Gvs,posclqpep)
cmpclq2=components(Gvsclq2)
Gvsclq2=induced_subgraph(Gvsclq2,names(cmpclq2$membership)[cmpclq2$membership==1])
vGvsclq2=names(V(Gvsclq2))
dGvsclq2=degree(Gvsclq2)
Gvclq2Louv=treeCrawLouv(vGvsclq2, Gvsclq2, s=10) #recursive Louvain clustering - yields a tree
Gvclq2Louv=lapply(Gvclq2Louv, function(x) {
  x=unlist(x)
  names(x)=NULL
  return(x)
})

rownames(Groupings)=pats
x=melt(Gvclq2Louv)
nx=x$value 
x=x$L1
names(x)=nx
rm(nx)  
co==c(11,2,)
plot(Gvsclq2, vertex.size=2.5, vertex.label=NA, vertex.color=x[vGvsclq2], palette=palette())


# Peptide map of the ClqClq graph clusters ------------------------
clq2ppsimple=aggregate(clqclqpp$Clq2, by=list(clqclqpp$Pep), "list")
x=clq2ppsimple$x
names(x)=clq2ppsimple$Group.1
clq2ppsimple=x
clq2ppsimple=sapply(clq2ppsimple,unique)

Clq2pp=aggregate(clqclqpp$Pep, by=list(clqclqpp$Clq2), unique)
x=Clq2pp$x
names(x)=Clq2pp$Group.1
Clq2pp=x[order(as.numeric(Clq2pp$Group.1))]
names(Clq2pp)=1:23
Clq2Ags=lapply(Clq2pp,function(l){
  l=unlist(l)
  X=chisq.posthoc.test(table(pepiused,names(pepiused) %in% l),simulate.p.value=T)
  n=nrow(X)
  i=1:(n/2)
  j=X[2*i,4]<0.05
  X=X[(2*i)-1,c(1,4)]
  X=X[order(X[,2], decreasing = T),]
  X[X[,2]>0,]
})
sink(file="Clq2pp_positiveAgsignif_new.txt")
print(Clq2Ags)
for (i in seq_along(Clq2pp)){
  l=Clq2pp[[i]]
  print(list(i,table(pepiused,names(pepiused) %in% l)))
}
sink()

Clq2Ags_m=melt(Clq2Ags)
Clq2Ags_m$Dimension=shortNames[Clq2Ags_m$Dimension]
Clq2Ags_m$Dimension[is.na(Clq2Ags_m$Dimension)]="Linker"
Clq2Ags_mx=acast(data=Clq2Ags_m, Dimension~L1)
Clq2Ags_mx[is.na(Clq2Ags_mx)]=0
heatmap3(Clq2Ags_mx, col=cpl1(1024),method="ward")
z=rowsums(Clq2Ags_mx)
names(z)=rownames(Clq2Ags_mx)
sink(file="clq2_Agsorted.txt")
print(sort(z, decreasing = T))
sink()
sink("Ags_excluded.txt")
print(pepi[!(shortNames[pepi] %in% names(z))])
sink()


plot(hclust((dist(t(Clq2Ags_mx)))))

X=Clq2Ags_mx
X[X==0]=rnorm(sum(X==0),0,1e-16)
X=X-min(X)+1e-17
X=log10(abs(X))
x=prcomp(X)
n=rownames(X) %in% viralshort
n[n]="Viral"
n[n!="Viral"]="Self"
fviz_screeplot(x, ncp=21)
fviz_pca_biplot(x, axes=c(1,2),repel=T, habillage=n)
fviz_pca_biplot(x, axes=c(3,4),repel=T, habillage=n)
fviz_pca_biplot(x, axes=c(5,6),repel=T, habillage=n)
fviz_pca_biplot(x, axes=c(7,8),repel=T, habillage=n)
sink(file="shortNm.txt")
print(shortNames)
sink()

Clq2Vir=sapply(Clq2pp, function(x){
  tb=c(0,0)
  names(tb)=c("Self","Vir")
  x=x %in% vir
  x[x]="Vir"
  x[!x=="Vir"]="Self"
  ti=table(x)
  tb[names(ti)]=ti
  return(tb)
})

chisq.posthoc.test(Clq2Vir, simulate.p.value=T)

clqVir=t(sapply(clqord, function(cq){
  cq=names(cq)
  l=length(cq)
  pos=sum(cq %in% vir)
  c(l-pos,pos)
}))

clqVir=aggregate(clqVir, by=list(lengths(clqord)), "sum")
chisq.posthoc.test(clqVir[,2:3], simulate.p.value=T)
xm=colSums(clqVir)
allmn=xm[3]/(sum(xm[2:3]))
x=clqVir$Group.1
y=clqVir$V2/(clqVir$V1+clqVir$V2)
plot(x, y, xlim=range(x), ylim=range(y), xlab="clique Size", ylab="Proportion Viral Epitopes") 
par(new=T)
plot(c(3,39),c(allmn,allmn), ty="l", xlim=range(x), ylim=range(y), xlab="", ylab="")


# MEME motifs of Clq2 peptides --------------------------------------------

clq2motifs=read.csv(file="XSTREMEclq2MOTIFS/motifs.csv")
x1=sapply(clq2ppsimple,function(z) z[1])
x2=sapply(clq2ppsimple,function(z) if (length(z)>1) return(z[2]) else return(NA))
x2=x2[!is.na(x2)]
clq2motifs=data.frame(clq2motifs,Clique=x1[clq2motifs$Peptides], Clique2=x2[clq2motifs$Peptides])
clq2motifs=clq2motifs[,-6]
Gvmeme=induced_subgraph(Gvs, clq2motifs$Peptides)
Gvmeme=set.vertex.attribute(Gvmeme, "motif", clq2motifs$Peptides, value=clq2motifs$Motif)
vGvmeme=names(V(Gvmeme))
################################################



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
flnm="clqFs.txt"

clqFs=lapply(1:21, function(pat){
  print(pat)
  proct=proc.time()
  clqGvs3_4=max_cliques(Gi[[pat]], min=3, max=4)
  J=cut(seq_along(clqGvs3_4), 45, labels=F)
  cl=makeCluster(15)
  clusterExport(cl, c("pat","mx","dgnf","N", "clucri","connfix","dunnfix","bhgamfix","BHgamma","clqGvs3_4","J"), envir = environment())
  diaclq=pbsapply(1:45, function(j){
    iJ=J==j
    bsmx=array(sample(mx[,-pat]), dim = dim(mx[,-pat]))
    sapply(clqGvs3_4[iJ], function(cq){
          cq=names(cq) 
          print(c(j,cq))
          n=length(cq)
          y=mx[cq,-pat]
          x=clucri(y, dgnf[-pat])
          bsx=sapply(1:500,function(i) {
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
  print(proc.time()-proct)
  nam=paste(pat,cls,"\n",sep = "_")
  write(nam, file=flnm, append=T)
  return(cls)
})

clqFsFromFile=read.table(file="clqFs.txt", header = F, row.names = NULL)
clqFs1=t(sapply(clqFsFromFile$V1,function(x) unlist(strsplit(x, split="_"))))
clqFs1=aggregate(clqFs1[,2], by=list(clqFs1[,1]), "list")$x

cofmxclqFs=SVMforValidation(clqFs1,X=mx,dgnf )
clqFTb=table(unlist(clqFs1))

y0=t(sapply(0:20,function(i){
  y=names(clqFTb)[clqFTb>i]
  if (length(y)>1) {
    plotMDS(mx[y,], col=dgnf, main=i)
    return(c(i,length(y),clucri(vs10[y,], dgnf)))
  }
  else {
    return(c(i,length(y),0))
  }
}))

y=names(clqFTb)[clqFTb>4]
plotMDS(mx[y,], col=dgnf, labels=NULL, pch=16)
callsclqf=names(clqFTb)[clqFTb>4]
#save(callsclqf, file="callsclqfs_nestedANOVA")

pdf(file="hmclqFall_51.pdf")
hm=heatmap.2(vs10[callsclqf,], hclustfun = hclwrd,  na.rm=F,key.title = NA,symkey=FALSE, cexCol = 0.7,cexRow = 0.5,
             trace="none",col=cpl1(1000),colsep=c(5,10,15),rowsep=c(5,10,15,28,36,40,47),
             margins = c(5,12), lwid=c(0.5,1), lhei = c(0.2,1))    
dev.off()
clcalls=callsclqf[rev(hm$rowInd)]
patord=hm$colInd
cli=cut(seq_along(clcalls),breaks = c(0,5,10,15,28,36,40,47,70), labels = F)
pati=cut(1:21,breaks = c(0,5,10,15,22), labels = F)

claov=sapply(1:8,function(i){
  n=length(cli[cli==i])
  y=c(mx[clcalls[cli==i],patord])
  cl=as.factor(pati)
  p=as.factor(rep(1:21,each=n))
  y=aggregate(y,by=list(p), "mean")$x
  a=aov(y~cl)
  x=summary(a)
  ta=TukeyHSD(a,"cl")
  ta=ta$cl[,1]*(ta$cl[,4]<0.01)
  return(ta)
})
closs=apply(claov,2,function(j) any(j[c(3,5,6)]>0))
clain=apply(claov,2,function(j) any(j[c(3,5,6)]<0))

x=clcalls
protclcls=sapply(x, function(p){
  pepiused[names(pepiused)==p]
})
tclcls=table(unlist(protclcls),cli)
chisq.test(tclcls, simulate.p.value = T)

aggainloss=c(-1,0,0,-1,1,-1,-1,-1)
names(aggainloss)=1:8
aggainloss=as.factor(aggainloss[cli])
tclglss=table(unlist(protclcls),aggainloss)
chisq.posthoc.test(tclglss, simulate.p.value = T)

virclcls=clcalls %in% vir
tvirclcls=table(virclcls,cli)
chisq.posthoc.test(tvirclcls, simulate.p.value=T)
x=clcalls[cli==5]
pepiused[x]

# Increased reactivity to viruses correlates with tumors - 
# EBV, papilloma and KSV, in the same cluster there are 2 epitopes from
# stromelysin and one each from claudin-6 and tyrosinase

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

clcallsmx=vs10[clcalls,]
Y=melt(clcallsmx)
Gr=apply(Groupings, 2, function(l) rep(l,each=51))
cls=rep(cli,21)
Y=data.frame(Y,Gr,Cluster=cls)
Y[,c(1,2,4:11)]=apply(Y[,c(1,2,4:11)],2,as.factor)
lmclcallsGr=lme(value~A*Cluster, random=~1|Var2, data=Y)
res=lsmeans(lmclcallsGr, pairwise~A|Cluster)
reres=summary(res$contrasts)
lmclcallsGr=lme(value~B*Cluster, random=~1|Var2, data=Y)
res=lsmeans(lmclcallsGr, pairwise~B|Cluster)
reres=summary(res$contrasts)
lmclcallsGr=lme(value~sex*Cluster, random=~1|Var2, data=Y)
res=lsmeans(lmclcallsGr, pairwise~sex|Cluster)
reres=summary(res$contrasts)

######### Comparisons -----------------------------------------------------
pp=names(pepiused)[names(pepiused) %in% pep]
calls=vs10[pp,]
x=pepiused[pp]
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

tGcsLC=table(GvcL,names(GvcL) %in% callsclqf)
chisqGvsLC=chisq.posthoc.test(tGcsLC, simulate.p.value = T)

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

q7cut=cut(q7ppp[names(GvcL)], c(-1,0.5,1.5,100), labels = F)
ti7GLC=table(GvcL,q7cut)

chisqi7GLC=chisq.posthoc.test(ti7GLC, simulate.p.value = T)
clqclqppi7=cbind(clqclqpp,Idiotypes=q7ppp[clqclqpp$Pep])
chisq.posthoc.test(table(clqclqppi7$Clq2, cut(clqclqppi7$Idiotypes, c(-1,0.5,1.5,100), labels = F)), simulate.p.value=T)

# Clusters of cliques 6,10 and 14 have a significant increase in idiotypes


# Correlation with viral

tvGLC=table(GvcL,vGvs %in% vir)

chisqvGLC=chisq.posthoc.test(tvGLC, simulate.p.value = T)

pepvircl6=vGvs[GvcL==6 & vGvs %in% vir]
length(unique(pepinfo$protein[pepinfo$peptide %in% pepvircl6]))
length(vnms)


# Correlations with NN Dia Groups -----------------------------------------



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

aovsex=sapply(1:14, function(i){
  print(i)
  m=vb10[Lvcl==i,]
  m=melt(m)
  m=cbind(m, Sex=sex[m$Var2], Pat=rep(pats,each=sum(Lvcl==i)))
  w=summary(aov(data=m, value~Sex+Error(Pat)))
  return(w[[1]][[1]][1,5])
})
aovsex=p.adjust(aovsex)


# Cor. Sex - NN - egonetwork ----------------------------------------------


# Gender specific reactivities correlate with antiviral in women
# and are anticorrelated in men


# Highly common idiotope determinants are crossreactive with both male and 
# female specific reactivties (yellow) but are not indifferent to gender 
# separataion. Unique idiotype reactivities are more common in man than in
# women.


# When using only the central vertex instead of the whole neighborhood
# these findings do not reach statistical significance


# ABO NN mean--------------------------------------------------------------

cl=makeCluster(15)
clusterExport(cl,c("vs10","BlGr","Gvs"))
clusterEvalQ(cl,{require(nlme)
                 require(igraph)
                 require(lsmeans())})
AB0mNN=t(pbsapply(vGvs, function(p){
    pp=names(ego(Gvs,1, p)[[1]])
    x=vs10[pp,]
    #print(which(vGvs==p))
    if (!is.null(nrow(x))) n=nrow(x) else n=1
    Y=c(x)
    pat=as.factor(rep(1:21,each=n))
    xA=as.factor(rep(BlGr$A,each=n))
    xB=as.factor(rep(BlGr$B,each=n))

    lmrA=lme(Y~xA,random=~1|pat)
    resA=lsmeans(lmrA, pairwise~xA)
    xA=summary(resA)
    lmrB=lme(Y~xB,random=~1|pat)
    resB=lsmeans(lmrB, pairwise~xB)
    xB=summary(resB)
    
    if(n>1) {y=data.frame(F_A=as.numeric(xA$contrasts[[2]]),
                 F_B=as.numeric(xB$contrasts[[2]]),
                 PadjA=as.numeric(xA$contrasts[[6]]),
                 PadjB=as.numeric(xA$contrasts[[6]]))
              return(y)
    }
    else {
      return(data.frame(F_A=0,F_B=0,PadjA=1,PadjB=1))
    }
   
}, cl=cl))
stopCluster(cl)

mx=vb10-min(vb10)+1

AB0mNN[,3:4]=apply(AB0mNN[,3:4],2,p.adjust,method="BH")

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
prof=as.factor(BlGr$A)
shamFs=sapply(1:50, function(i){
  print(i)
  if (i>1) prof=sample(as.factor(BlGr$A))
  crit=sapply(vGvs, function(p){
    pp=names(ego(Gvs,1,p)[[1]])
    x=mx[pp,]
    if (!is.null(nrow(x))) {
      n=nrow(x)
      Y=colMeans(x)
    }
    else {
      n=1
      Y=x
    }
    xA=as.factor(prof)
    diff(aggregate(Y, by=list(xA),"mean")[,2])
  })
  return(crit)
})

nmatr=paste("X",i,sep="",collapse="")
Gvs=set.vertex.attribute(Gvs, nmatr, value = crit)
write.graph(Gvs, format = "graphml", file="Gvs.graphml")

for (i in 1:10){
  nmatr=paste("X",i,sep="",collapse="")
  Gvs=delete_vertex_attr(Gvs, nmatr)
}

plot(AB0loNN[,2]~AB0loNN[,1])
plot(lm(AB0loNN[,3]~AB0loNN[,1])$residuals~AB0loNN[,2])

nnq10=table(cut(AB0loNN[,1],c(-10,-0.152,0.152,10),labels = c("Aneg","Anone","Apos")),
            cut(AB0loNN[,2],c(-10,-0.152,0.152,10),labels = c("Bneg","Bnone","Bpos")))
chisq.posthoc.test(nnq10, simulate.p.value=T)
chisq.test(nnq10, simulate.p.value=T)

Zass_pro(prof=BlGr$A)


Sdist=t(sapply(vGvs, function(p){
  pp=names(ego(Gvs,1, p)[[1]])
  x=vs10[pp,]
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

# PCA cor

x=prcomp(vs10)
Groupings=cbind(BlGr,dgnC,dgnG,dgnM)
PCAcor=sapply(1:5,function(i) {
  print(i)
  y=sapply(1:6,function(j){
    w=wilcox.test(x$rotation[,i]~Groupings[,j])
    ef=diff(aggregate(x$rotation[,i], by=list(Groupings[,j]), "mean")$x)
    p=w$p.value
    if (p>0.1) p=1 
    return(c(p,ef))
  })
  y=data.frame(t(y[1,]),t(y[2,]))
  cnm=c("A","B","Sex","Contr","GBM","Meta")
  cnm=c(cnm,paste("eff",cnm,sep="_"))
  colnames(y)=cnm
  return(y)
})

PCAcor[1:6,]=p.adjust(PCAcor[1:6,])


# New ABO Louvain ---------------------------------------------------------
n=length(GvcL)
AB0xLouv=melt(vs10[names(GvcL),])
grs=apply(Groupings,2,function(l) rep(l,each=n))
grs=as.data.frame(grs)
Grrep=grs
Grrep$A=gsub(TRUE, "A", Grrep$A)
Grrep$A=gsub(FALSE, "nonA", Grrep$A)
Grrep$B=gsub(TRUE, "B", Grrep$B)
Grrep$B=gsub(FALSE, "nonB", Grrep$B)
Grrep$N=gsub(TRUE, "O", Grrep$N)
Grrep$N=gsub(FALSE, "nonO", Grrep$N)
Grrep$dgnC=gsub(TRUE, "Contr", Grrep$dgnC)
Grrep$dgnC=gsub(FALSE, "Tu", Grrep$dgnC)
Grrep$dgnG=gsub(TRUE, "GBM", Grrep$dgnG)
Grrep$dgnG=gsub(FALSE, "nonGBM", Grrep$dgnG)
Grrep$dgnM=gsub(TRUE, "Meta", Grrep$dgnM)
Grrep$dgnM=gsub(FALSE, "nonMeta", Grrep$dgnM)
Grrep$sex=gsub(TRUE, "F", Grrep$sex)
Grrep$sex=gsub(FALSE, "M", Grrep$sex)
AB0xLouv=data.frame(Pat=as.factor(AB0xLouv$Var2), Value=AB0xLouv$value, 
             Louv=as.factor(rep(GvcL,21)),A=as.factor(Grrep$A), B=as.factor(Grrep$B), 
             O=as.factor(Grrep$N),Contr=as.factor(Grrep$dgnC), 
             GBM=as.factor(Grrep$dgnG),Meta=as.factor(Grrep$dgnM),Sex=as.factor(Grrep$sex))
ABLlm=lme(data=AB0xLouv, Value~A*Sex*Louv,random=~1|Pat)
resABLlm=lsmeans(ABLlm, pairwise~A|Sex|Louv)
resABLlm=summary(resABLlm)
sink(file="blGrA_SeX_Louv.txt")
print(resABLlm)
sink()
ABLlm=lme(data=AB0xLouv, Value~B*Sex*Louv,random=~1|Pat)
resABLlm=lsmeans(ABLlm, pairwise~B|Sex|Louv)
resABLlm=summary(resABLlm)
sink(file="blGrB_SeX_Louv.txt")
print(resABLlm)
sink()

DiaCLlm=lme(data=AB0xLouv, Value~Contr*Louv,random=~1|Pat)
resDiaLlm=lsmeans(DiaCLlm, pairwise~Contr|Louv)
resDiaLlm=summary(resDiaLlm)
reslmC=t(sapply(1:14, function(i) resDiaLlm$contrasts[i,]))
reslmC=reslmC[order(unlist(reslmC[,3])),]
write.csv(reslmC, file="reslmC.csv")

DiaGLlm=lme(data=AB0xLouv, Value~GBM*Louv,random=~1|Pat)
resDiaLlm=lsmeans(DiaGLlm, pairwise~GBM|Louv)
resDiaLlm=summary(resDiaLlm)
reslmG=t(sapply(1:14, function(i) resDiaLlm$contrasts[i,]))
reslmG=reslmG[order(unlist(reslmG[,3])),]
write.csv(reslmG, file="reslmG.csv")

DiaMLlm=lme(data=AB0xLouv, Value~Meta*Louv,random=~1|Pat)
resDiaLlm=lsmeans(DiaMLlm, pairwise~Meta|Louv)
resDiaLlm=summary(resDiaLlm)


# Modularities ------------------------------------------------------------

modtb=t(sapply(pepi,function(p){
  Zmod(Gvs, as.numeric(vGvs %in% names(pepiused)[pepiused==p])+1)
}))
colnames(modtb)=c("Modularity","Z-score")
modv=Zmod(Gvs, as.numeric(vGvs %in% vir)+1)
modlouv=Zmod(Gvs, as.numeric(vGvs %in% vir)+1)
modtb=rbind(modtb,`All Viral`=modv)
x=GvcL[vGvs]
x[is.na(x)]=100
modlouv=Zmod(Gvs, as.numeric(x))
modtb=rbind(modtb,`Louvain Clustering`=modlouv)
modtb=modtb[order(modtb[,2], decreasing = T),]
modtb=cbind(modtb,Viral=rownames(modtb) %in% pepi[viral])
x=sapply(rownames(modtb),function(x) {
  if (x %in% pepi) {
    mean(dGvs[names(pepiused)[pepiused %in% x]])
  }
  else {
    return(0)
  }
})
x[1]=mean(dGvs)
xv=dGvs[vir[!is.na(vir)]]
x[6]=mean(xv[!is.na(xv)])
modtb=cbind(modtb,`Mean Degree`=x)

write.table(modtb,file = "modtb.csv")
### sex


### degree
degtr=cut(dGvs, quantile(dGvs, c(0,0.33,0.67,1))-c(0.001,0,0,-0.001), labels=F)
degass=Zass(Gvs, degtr)


# AB0

