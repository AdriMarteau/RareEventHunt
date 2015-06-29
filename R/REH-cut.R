tree.fRate <- function(set.sub, set.full)
{
  if(length(set.sub)==0) return(0)
  if(length(intersect(set.full,set.sub))==0) return(0)
  length(intersect(intersect(set.full,set.sub),frdUsed))/length(intersect(set.full,set.sub))
}

tree.ProbaAtSplit <- function(Data, Cut, set.full)
{
  sub.set=which(Data<=Cut)
  p1=tree.fRate(sub.set, set.full)
  p2=tree.fRate(setdiff(set.full,sub.set), set.full)
  n1=length(intersect(set.full,sub.set))
  n2=length(intersect(set.full,setdiff(set.full,sub.set)))
  score=p1-p2
  c(p1,p2,score,n1,n2)
}

tree.UncertRedCut <- function(Data, Reg.cuts, Reg.names, set.full)
{
  allcutres = lapply(Reg.names, 
                     function(s){
                       sapply(Reg.cuts[[s]], function(x) 
                       {
                         x=tree.ProbaAtSplit(Data[[s]],x,set.full); 
                         if(x[4]*x[5]==0 || min(x[1:2])<prb.lim.sensi || max(x[1:2])>1-prb.lim.sensi) return(0)
                         return(sign(x[3])/(1+max(x[1]*(1-x[1])/x[4],x[2]*(1-x[2])/x[5])))
                       })
                     }
  )
  names(allcutres)<-Reg.names
  listbest = sapply(Reg.names, function(s){ u=which.max(abs(allcutres[[s]])); return(c(u,allcutres[[s]][u]))})
  bestvar = which.max(abs(listbest[2,]))
  x=listbest[[1,bestvar]]
  sgn = sign(listbest[[2,bestvar]])
  best.reg = names(bestvar)
  best.cut = Reg.cuts[[best.reg]][x]
  p=tree.ProbaAtSplit(Data[[best.reg]],best.cut,set.full)
  c(best.reg,best.cut,sgn,p[1],p[2],p[4],p[5])
}

tree.ObjectiveCatExists <- function(Data, Reg.cuts, Reg.names, set.full)
{
  allcutres = lapply(Reg.names, 
                     function(s){
                       sapply(Reg.cuts[[s]], function(x) 
                       {
                         x=tree.ProbaAtSplit(Data[[s]],x,set.full); 
                         if(x[3]>0) return((x[1]>obj.min.prb && x[4] > obj.min.size)*x[1]*x[4])
                         return(-(x[2]>obj.min.prb && x[5] > obj.min.size)*x[2]*x[5])
                       })
                     }
  )
  names(allcutres)<-Reg.names
  listbest = sapply(Reg.names, function(s){ u=which.max(abs(allcutres[[s]])); return(c(u,allcutres[[s]][u]))})
  if(max(abs(listbest[2,]))==0) return(NULL)
  bestvar = which.max(abs(listbest[2,]))
  x=listbest[[1,bestvar]]
  sgn = sign(listbest[[2,bestvar]])
  best.reg = names(bestvar)
  best.cut = Reg.cuts[[best.reg]][x]
  p=tree.ProbaAtSplit(Data[[best.reg]],best.cut,set.full)
  c(best.reg,best.cut,sgn,p[1],p[2],p[4],p[5])
}

tree.LargeSampleCut <- function(Data, Reg.cuts, Reg.names, set.full)
{
  allcutres = lapply(Reg.names, 
                     function(s){
                       sapply(Reg.cuts[[s]], function(x) 
                       {
                         x=tree.ProbaAtSplit(Data[[s]],x,set.full); 
                         if(x[4]*x[5]==0 || min(x[4],x[5]) < rho_*(x[4]+x[5])) return(0)
                         if(x[3] >0) return(x[3]*sqrt(x[4]))
                         else  return(x[3]*sqrt(x[5]))
                       })
                     }
  )
  names(allcutres)<-Reg.names
  listbest = sapply(Reg.names, function(s){ u=which.max(abs(allcutres[[s]])); return(c(u,allcutres[[s]][u]))})
  bestvar = which.max(abs(listbest[2,]))
  x=listbest[[1,bestvar]]
  sgn = sign(listbest[[2,bestvar]])
  best.reg = names(bestvar)
  best.cut = Reg.cuts[[best.reg]][x]
  p=tree.ProbaAtSplit(Data[[best.reg]],best.cut,set.full)
  c(best.reg,best.cut,sgn,p[1],p[2],p[4],p[5])
}

tree.SmallSampleCut <- function(Data, Reg.cuts, Reg.names, set.full)
{
  allcutres = lapply(Reg.names, 
                     function(s){
                       sapply(Reg.cuts[[s]], function(x) 
                       {
                         x=tree.ProbaAtSplit(Data[[s]],x,set.full); 
                         if(x[4]*x[5]==0) return(0)
                         if(x[3] >0) return(x[3]*sqrt(x[4]))
                         else  return(x[3]*sqrt(x[5]))
                       })
                     }
  )
  names(allcutres)<-Reg.names
  listbest = sapply(Reg.names, function(s){ u=which.max(abs(allcutres[[s]])); return(c(u,allcutres[[s]][u]))})
  bestvar = which.max(abs(listbest[2,]))
  x=listbest[[1,bestvar]]
  sgn = sign(listbest[[2,bestvar]])
  best.reg = names(bestvar)
  best.cut = Reg.cuts[[best.reg]][x]
  p=tree.ProbaAtSplit(Data[[best.reg]],best.cut,set.full)
  c(best.reg,best.cut,sgn,p[1],p[2],p[4],p[5])
}

tree.RefineCut <- function(Data, Reg.cuts, Reg.names, set.full)
{
  allcutres = sapply(Reg.names, 
                     function(s){
                       sapply(Reg.cuts[[s]], function(x) {tree.ProbaAtSplit(Data[[s]],x,set.full)})})
  listbest = sapply(Reg.names, function(s){ 
    if(!is.null(dim(allcutres[[s]][1:2,]))) x=apply(allcutres[[s]][1:2,],2,min)
    else x=min(allcutres[[s]][1:2,])
    y=which(allcutres[[s]][4,]*allcutres[[s]][5,]==0)
    x[y]<-Inf
    u=which.min(x)
    return(c(u,x[u]))}
  )	
  bestvar = which.min(abs(listbest[2,]))
  x=listbest[[1,bestvar]]
  sgn = sign(listbest[[2,bestvar]])
  best.reg = names(bestvar)
  best.cut = Reg.cuts[[best.reg]][x]
  p=tree.ProbaAtSplit(Data[[best.reg]],best.cut,set.full)
  c(best.reg,best.cut,sgn,p[1],p[2],p[4],p[5])
}