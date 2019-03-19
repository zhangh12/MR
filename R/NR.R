
NR <- function(para, map, the0, inv.the, gam0, inv.gam, 
               n1, n0, ref){
  
  f <- score.gim(para, map, the0, inv.the, gam0, inv.gam, 
                 n1, n0, ref)
  para - f
  
}
