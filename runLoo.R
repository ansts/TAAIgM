i=commandArgs(trailingOnly = TRUE)
i=as.double(i)

source("bhgamfix.R")
source("BHgamma.R")
source("connfix.R")
source("dunnfix.R")
source("rfe.R")
source("Loo.R")
load("vars")

slf=Loo(i)

save(slf, file=paste("Fs","part",i,sep="_",collapse=""))
