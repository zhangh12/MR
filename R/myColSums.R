
## Wheeler Bill developed this code, I did not change it

myColSums <- function(x, nr, nc, na.rm=FALSE) {
  z <- .Internal(colSums(x, nr, nc, na.rm))  
  z
} # END: myColSums

