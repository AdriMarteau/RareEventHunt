
printprb <- function(p) {
  s=paste("     ",format(p*100, digits = 2, nsmall = 2),sep="")
  substr(s,nchar(s)-5,nchar(s))
}

getSplitInfo <- function(myTree, info, leaf.ret)
{
  if(attr(myTree, "members")==1) return(leaf.ret)
  return(c(attr(myTree, info),Recall(myTree[[1]],info,leaf.ret),Recall(myTree[[2]],info,leaf.ret)))
}

getNodeType <- function(myTree,tp="N")
{
  if(attr(myTree, "members")==1) return(tp)
  tp="N";	return(c(tp,Recall(myTree[[1]],"L"),Recall(myTree[[2]],"R")))
}

getLeafInfo <- function(myTree, info, val)
{
  if(attr(myTree, "members")==1) return(val)
  l.val=attr(myTree, paste(info,"l",sep="."))
  r.val=attr(myTree, paste(info,"r",sep="."))
  return(c(NA,Recall(myTree[[1]],info,l.val),Recall(myTree[[2]],info,r.val)))
}

runRandom <- function(q__, plot.res=TRUE)
{
  q_ <<- q__
  Debug<<-0
  k <<- sample(Reg.names,q_,replace=F)
  myTree<<-tree.built(Data, Reg.cuts, k, set.full, sensi.a, sensi.p, sensi.q, min.p.l,maxDepth=maxDepth)
  if(plot.res) g.plot(myTree, main=paste(q_[1], "params:", paste(k, collapse=", ")),conf.r=3)
}

RerunRandom <- function(...)
{
  myTree<<-tree.built(Data, Reg.cuts, k, set.full, sensi.a, sensi.p, sensi.q, min.p.l,maxDepth=maxDepth)
  g.plot(myTree, main=paste(q_[1], "params:", paste(k, collapse=", ")),conf.r=3)
}


cor.cat.cut=c("black","red","green","blue")


g.plot <- function(myTree,horizontal=T,main="",conf.r=1)
{	
  cuts.lab=getSplitInfo(myTree,"split.label", "")
  cuts.me=getSplitInfo(myTree,"members",0)
  cuts.type = getSplitInfo(myTree,"cut.type",0)
  cat.prb = getLeafInfo(myTree,"p",0)/100
  cat.n   = getLeafInfo(myTree,"n",0)
  
  cat.r = 1.96*sqrt(cat.prb*(1-cat.prb)/cat.n)*100
  cat.r = conf.r*cat.r/mean(cat.r,na.rm=T)/0.66
  offs =  1+cumsum(as.numeric(cuts.me==0))
  if(horizontal) 
  {
    op <- par(mar=c(1,1,1,1))
    side=4
    cuts.y=offs+getSplitInfo(myTree,"midpoint",-1)
    ylim=c(0,max(offs))
    cuts.x=getSplitInfo(myTree,"height",0)
    xlim=c(max(cuts.x)*(1.2),-0.5)
  } else {
    side=3
    op <- par(mar=c(1,1,1,1))
    cuts.x=offs+getSplitInfo(myTree,"midpoint",-1)
    xlim=c(0,max(offs))
    cuts.y=getSplitInfo(myTree,"height",0)
    ylim=c(-0.5,max(cuts.y)*(1.2))
  }
  plot(myTree, hor=horizontal, xlim=xlim, ylim=ylim, axes=F,main=main)
  points(cuts.x,cuts.y,col=cor.cat.cut[cuts.type+1],pch=16)
  points(cuts.x,cuts.y,cex=cat.r)
  text(cuts.x,cuts.y,cuts.lab,font=2,col="blue",pos=side)
  
  par(op)
} 
