
#' @title Calculate the sample skewness.
#' @description This function is used to compute the skewness coefficient of a sample.
#' @param x a numeric vector
#' @return a numeric vector of size 1
#' @examples
#' \dontrun{
#' x<-rnorm(20)
#' sk(x)
#' }
#' @export
sk <- function(x){
  mean <- mean(x)
  m3 <- mean((x-mean)^3)
  m2 <- mean((x-mean)^2)
  m3/m2^1.5
}

