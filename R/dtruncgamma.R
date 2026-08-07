#' truncated gamma density
#'
#' @param x point at which density is computed
#' @param shape shape parameter
#' @param scale scale parameter
#' @param a lower limit (default is 0)
#' @param b upper limit (default is Inf)
#'
#' @returns a real number
#' @export
#'
#' @examples
#' dtruncgamma(x=1,shape=2, scale=3, a=0, b=20)
dtruncgamma <- function(x, shape, scale, a = 0, b = Inf) {
  # Calculate the standard Gamma density
  dens <- stats::dgamma(x, shape = shape, scale = scale)

  # Calculate the normalization constant: Area under the curve between a and b
  norm_const <- stats::pgamma(b, shape = shape, scale = scale) - stats::pgamma(a, shape = shape, scale = scale)

  # Normalize and set values outside the interval [a, b] to 0
  result <- dens / norm_const
  result[x < a | x > b] <- 0

  return(result)
}
