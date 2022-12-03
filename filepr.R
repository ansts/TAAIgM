#
# filepr(ocessing) recovers all files from path and applies grpLocNor 
# to the in parallel batch (if sh(ow) is F) or one by one if sh=T to show 
# the graphic presentations
#

filepr<-function(path, sh=F, RG=T) { 
    require(parallel)
    lf=list.files(path = path)
    if (sh==F) {
      cl <- makeCluster(4)
      ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
      clusterExport(cl, ex)
      parLapply(cl,lf, function(x){gprLocNorSVM(x, p=path, shw = sh, Rsw=RG)})
      stopCluster(cl)
    }
    else {
      lapply(lf, function(x){gprLocNorSVM(x, p=path, shw=sh, Rsw=RG)})
    }
}