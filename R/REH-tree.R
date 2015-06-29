#source(file.path(REHlib.path,"REH.lib","REH_global.r"))
#source(file.path(REHlib.path,"REH.lib","REH_annex.r"))
#source(file.path(REHlib.path,"REH.lib","REH_scoring.r"))
#source(file.path(REHlib.path,"REH.lib","REH_plot.r"))


# +----------------------------+
# | Function to build the tree |
# +----------------------------+
REH <- function(Data, Reg.cuts, Reg.names, set.full, sensi.a, sensi.p, sensi.q, min.p.l=1, maxDepth=Inf)
{
  tree.built.env.tmp<-new.env(parent=globalenv())
  assign("r.leaf.sum", 0, envir=tree.built.env.tmp)
  assign("my.log", "", envir=tree.built.env.tmp)
  ret = tree.built.rec(Data=Data, Reg.cuts=Reg.cuts, Reg.names=Reg.names, set.full=set.full,
                       sensi.a=sensi.a, sensi.p=sensi.p, sensi.q=sensi.q,
                       min.p.l, maxDepth=maxDepth, 
                       cut.type=1, env=tree.built.env.tmp)
  RET = list(tree=ret,
             log  = get("my.log", envir=tree.built.env.tmp),
             eval = function(x) { TreeParse(x,ret)},
             xprt = function()  {TreeExport(ret,0)})
  return(RET)
}

TreeParse <- function(x, tree) 
{
  if(attr(tree,"height")==0) return(NA)
  var_ = attr(tree, "var")
  cut_ = attr(tree, "x.cut")
  sgn_ = attr(tree, "dir") > 0
  
  if((x[[var_]] <= cut_) == sgn_ & !is.na(x[[var_]]))
  {
    tmp<-Recall(x,tree[[2]])
    if(is.na(tmp)) return(attr(tree, "p.r")) else return(tmp)
  } else {
    tmp<-Recall(x,tree[[1]])
    if(is.na(tmp)) return(attr(tree, "p.l")) else return(tmp)
  }
}

TreeExport <- function(tree, i)
{
  if(attr(tree,"height")==0) return(NA)
  var_ = attr(tree, "var")
  cut_ = attr(tree, "x.cut")
  sgn_ = floor((as.numeric(attr(tree, "dir")) +1)/2)
  
  txt = paste(var_,cut_,sgn_,sep=":")
  txt = paste("<",txt,sep="|")
  
  i=1
  tmp.r<-Recall(tree[[2]],i)
  if(is.na(tmp.r)) txt = paste(txt,-1,sep="|") else { txt = paste(txt,i,sep="|");i=i+tmp.r$nline}
  
  tmp.l<-Recall(tree[[1]],i)
  if(is.na(tmp.l)) txt = paste(txt,-1,sep="|") else { txt = paste(txt,i,sep="|");i=i+tmp.l$nline}
  
  txt = paste(txt,attr(tree, "p.r"),attr(tree, "p.l"),">",sep="|")
  if(!is.na(tmp.r)) txt = paste(txt,tmp.r$txt,sep="\n")
  if(!is.na(tmp.l)) txt = paste(txt,tmp.l$txt,sep="\n")
  
  return(list(nline=i,txt =txt))
  
}



# +--------------------------------------+
# | Recursive function to build the tree |
# +--------------------------------------+
tree.built.rec <- function(Data, Reg.cuts, Reg.names, set.full,
                           sensi.a, sensi.p, sensi.q, 
                           min.p.l, maxDepth, cur.depth = 1,
                           cut.type = 1, is.right=F, env=env)
{
  Debug<<-set.full
  if(maxDepth==-1) return(NULL) #stop
  x = tree.ObjectiveCatExists(Data, Reg.cuts, Reg.names, set.full)
  n=length(set.full)
  N=dim(Data)[1]
  if(is.null(x)) {
    if(cut.type==1 & n<sensi.q*N) cut.type=2
    x=tree.OptimalCut(Data, Reg.cuts, Reg.names, set.full, cut.type)
  } else {
    cut.type = 3
  }
  cTmp=as.numeric(x[2])
  sgn=x[3]
  if(x[3]>0)
  {
    p.r=as.numeric(x[4]);	n.r=as.numeric(x[6]);	p.l=as.numeric(x[5]);	n.l=as.numeric(x[7])
    sgn.str = "<="
    sub.set.tmp=intersect(set.full,which(Data[[x[1]]]<=cTmp))
  } else{
    p.l=as.numeric(x[4]);	n.l=as.numeric(x[6]);	p.r=as.numeric(x[5]);	n.r=as.numeric(x[7])
    sgn.str = " >"
    sub.set.tmp=intersect(set.full,which(Data[[x[1]]]>cTmp))
  }
  n=length(sub.set.tmp) + get("r.leaf.sum", envir=env)
  r.child = NULL;	l.child = NULL
  my.log.tmp = get("my.log", envir=env)
  my.log.tmp = paste(my.log.tmp,cut.type.string[cut.type],sep="")
  my.log.tmp = paste(my.log.tmp,"cut:",	substr(paste(x[1],"     ",sep=""),1,5), sgn.str,
                     substr(paste(x[2],"     ",sep=""),1,5), 
                     "[p.r=", printprb(p.r), "| p.l=", printprb(p.l),	"]\n")
  
  if(length(sub.set.tmp)==length(set.full) | length(sub.set.tmp)==0) return(NULL)
  assign("my.log", my.log.tmp, envir=env)
  if(n >= N*(sensi.a + sensi.p) && cut.type!=3) {	
    r.child = Recall(Data, Reg.cuts, Reg.names, sub.set.tmp ,
                     sensi.a, sensi.p, sensi.q, 
                     min.p.l, maxDepth-1, cur.depth = cur.depth + 1,
                     is.right=T, env=env)
    if(p.l > min.p.l) { 
      #will refine even if the max deapth if reached...	
      sub.set.tmp2 = setdiff(set.full,sub.set.tmp)
      l.child = Recall(Data, Reg.cuts, Reg.names, sub.set.tmp2,
                       sensi.a, sensi.p, sensi.q, 
                       min.p.l, maxDepth=0, cur.depth = cur.depth + 1,
                       is.right=F, cut.type = 1, env=env)
    }	
  } 
  if(n <  N*(sensi.a) && cut.type!=3) {	
    sub.set.tmp2 = setdiff(set.full,sub.set.tmp)
    l.child = Recall(Data, Reg.cuts, Reg.names, sub.set.tmp2,
                     sensi.a, sensi.p, sensi.q, 
                     min.p.l, maxDepth-1, cur.depth = cur.depth + 1,
                     is.right=F, env=env)
  }
  
  p.r=p.r*100; p.l=p.l*100
  ddtmp = list()
  ddtmp$merge=matrix(c(-1, -2), nc=2, byrow=TRUE )
  ddtmp$height <- 1
  ddtmp$labels <- c(paste(format(p.l, digits = 2, nsmall = 2),"%", " (",n.l,")",sep=""),
                    paste(format(p.r, digits = 2, nsmall = 2),"%", " (",n.r,")",sep=""))
  ddtmp$order <- 1:2
  class(ddtmp) <- "hclust"
  ddtmp <- as.dendrogram(ddtmp)
  
  if(length(r.child)==0 & length(l.child)==0)
  {
    if(is.right) 
    {
      assign("r.leaf.sum"   , n, envir=env)
    }
  } else {
    if(length(r.child)==0)	ddtmp= merge(l.child,ddtmp[[2]])
    else if(length(l.child)==0)	ddtmp= merge(ddtmp[[1]],r.child)
    else ddtmp= merge(l.child,r.child)
  }	
  
  attr(ddtmp,"split.label") <- paste(x[1],sgn.str,x[2])
  attr(ddtmp,"var") <- x[1]
  attr(ddtmp,"dir") <- sgn
  attr(ddtmp,"x.cut") <- cTmp
  attr(ddtmp,"p.r") <- p.r
  attr(ddtmp,"n.r") <- n.r
  attr(ddtmp,"p.l") <- p.l
  attr(ddtmp,"n.l") <- n.l
  attr(ddtmp,"cut.type") <- cut.type
  
  return(ddtmp)	
}