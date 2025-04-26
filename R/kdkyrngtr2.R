#' compute max(kD) and max(kY) for total R2-based analysis
#'
#' @param data data frame for analysis
#' @param outcome name of outcome variable
#' @param treatment name of treatment variable
#' @param bnch_reg name(s) of benchmark covariate(s)
#' @param other_reg name(s) of other covariates
#'
#' @returns a data frame with 2 columns and 1 row:
#' \item{kd_high}{max(kD), a scalar}
#' \item{ky_high}{max(kY), a scalar}
#' @export
#'
#' @examples
#' require("sensemakr")
#' Y <- "peacefactor"
#' D <- "directlyharmed"
#' X <- "female"
#' X_oth <- c("village","age","farmer_dar","herder_dar","pastvoted","hhsize_darfur")
#'
#'
#' r1 <- kdkyrngtr2(data=darfur,outcome=Y,treatment=D,bnch_reg=X,other_reg=X_oth)
#'
kdkyrngtr2 <- function(
    data,outcome,treatment,bnch_reg,other_reg=NULL
    ){

  # ------ Create data set
  d1 <- as.data.frame(data)

  # ------ Initial computations
  # --- R2DXj
  DXj <- stats::as.formula(paste(treatment, paste(bnch_reg, collapse= "+"), sep = "~"))
  R2DXj <- summary(stats::lm(DXj, data = d1))$r.squared

  # ---- R2DX
  DX <- stats::as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  R2DX <- summary(stats::lm(DX,data = d1))$r.squared

  # ---- R2YXj
  YXj <- stats::as.formula(paste(outcome, paste(bnch_reg, collapse= "+"), sep = "~"))
  R2YXj <- summary(stats::lm(YXj, data = d1))$r.squared

  # ----- R2YX
  YX <- stats::as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  R2YX <- summary(stats::lm(YX,data = d1))$r.squared


  # Return values
  return(
    data.frame(
      # Max of kD
      kd_high=((1-R2DX)/R2DXj),
      # Max of kY
      ky_high=((1-R2YX)/R2YXj)
    )
  )

}
