

inv.matrix <- function(mat){
  
  MinEig <- 0.001
  dec <- eigen(mat, symmetric=T)
  iv <- 1/(dec$values * (dec$values >= MinEig) + MinEig * (dec$values < MinEig))
  inv <- dec$vectors %*% diag(iv) %*% t(dec$vectors)
  
  inv
  
}
