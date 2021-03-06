#' @title  Random points from BS density.
#' @description  This function is used to generate random samples from Bart Simpson density.
#' @param n sample size
#' @return random sample from BS destiny
#' @examples
#' \dontrun{
#' n<-1000
#' y<-rBS(1000)
#' }
#' @export
rBS <- function(n){ 
  u <- runif(n)
  y <- u
  ind <- which(u > 0.5) #index for those generated from N(0,1)
  y[ind] <- rnorm(length(ind), 0, 1)
  for (j in 0:4) {
    ind <- which(u > j * 0.1 & u <= (j+1) * 0.1)
    #index for those generated from N(j/2-1,1/10^2)
    y[ind] <- rnorm(length(ind), j/2 -1, 1/10)
  }
  return(y)
}



