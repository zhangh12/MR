
hess.cc <- function(para, para.id, ref, the0, gam0, inv.V, nsample){

  h <- numDeriv::jacobian(score.cc, para, para.id = para.id, ref = ref, the0 = the0, gam0 = gam0, inv.V = inv.V, nsample = nsample)
  colnames(h) <- names(para)
  rownames(h) <- names(para)
  h

}
