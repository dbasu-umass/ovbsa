#' quasi-triangular probability distribution function
#'
#'`linvx()` is deprecated; please use [baci()] instead.
#'
#' @param x point (scalar) at which pdf is evaluated
#' @param xvec vector of all possible x values
#' @param k mode and median of the distribution
#'
#' @returns the value (scalar) of the pdf at x
#' @export
#'
#'@keywords internal
#'
#' @examples
#' xfull <- runif(n=100,min=0,max=10)
#' xpoint <- 5
#' xmod <- 2
#' res_pdf <- linvx(x=xpoint,xvec=xfull,k=xmod)
#'
linvx <- function(x,xvec,k){

  .Deprecated("baci", msg = "'linvx' is deprecated. Please use 'baci()' instead.")

  # right end point
  b <- max(xvec)

  # probability distribution function
  fx <- ifelse(
    x<=k,
    (1/k^2)*x,
    1/(2*log(b/k)*x)
  )

  # Return
  return(fx)

}
