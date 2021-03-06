#' @title Probability estimator of BS density.
#' @description This function is used to compute the estimated Bart Simpson density .
#' @param x estimate points
#' @param y samples
#' @param h bandwidth
#' @return estimated density at x
#' @examples
#' \dontrun{
#' x<-c(1,2,3)
#' y<-rnorm(20)
#' h<-0.5
#' peBS(x,y,h)
#' }
#' @export
peBS <- function(x, y, h){ # probability estimator of BS density

  m <- length(x)
  n <- length(y)
  ye <- rep(0, m)
  for (i in 1:n){
    ye <- ye + as.numeric((x >= y[i] - h) & (x < y[i] + h))
  }
  ye <- ye / (2*h*n)
  return(ye)
}


