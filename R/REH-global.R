prb.lim.sensi <<- 0.01/100
obj.min.prb   <<- 0.15
obj.min.size  <<- 150
rho_          <<- 0.1

cut.type.string = c("Large","Small","Optim") 

tree.OptimalCut <- function(Data, Reg.cuts, Reg.names, set.full, cut.type)
{
  if(cut.type==1) return(tree.LargeSampleCut(Data, Reg.cuts, Reg.names, set.full))
  if(cut.type==2) return(tree.UncertRedCut(Data, Reg.cuts, Reg.names, set.full))
  if(cut.type==-1) return(tree.RefineCut(Data, Reg.cuts, Reg.names, set.full))
}