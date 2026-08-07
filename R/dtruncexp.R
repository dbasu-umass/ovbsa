#' truncated exponential density
#'
#' @param x point between a and b at which density is computed
#' @param rate rate parameter
#' @param a lower bound
#' @param b upper bound
#'
#' @returns a real number
#' @importFrom stats dexp pexp
#' @export
#'
#' @examples
#' dtruncexp(x=1, rate=1, a=0, b=10)
dtruncexp <- function(x, rate = 1, a = 0, b = Inf) {
  # Calculate the standard exponential density
  dens <- stats::dexp(x, rate = rate)

  # Calculate the normalization constant: Area under the curve between a and b
  norm_const <- stats::pexp(b, rate = rate) - stats::pexp(a, rate = rate)

  # Normalize and set values outside the interval [a, b] to 0
  result <- dens / norm_const
  result[x < a | x > b] <- 0

  return(result)
}
