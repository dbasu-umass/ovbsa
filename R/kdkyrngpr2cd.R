#' compute max(kD) and max(kY) for partial R2-based analysis with conditioning on treatment
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
#' r1 <- kdkyrngpr2cd(data=darfur,outcome=Y,treatment=D,bnch_reg=X,other_reg=X_oth)
#'
kdkyrngpr2cd <- function(
    data,outcome,treatment,bnch_reg,
    other_reg=NULL){

  # ------ Create data set
  d1 <- as.data.frame(data)

  # ------ Initial computations

  # --- R2DXj.Xmj
  reg1 <- stats::as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  r1_d <- summary(stats::lm(reg1, data=d1))$r.squared

  if(is.null(other_reg)){
    r2_d <- 0
  } else{
    reg2 <- stats::as.formula(
      paste(treatment, paste(other_reg, collapse= "+"), sep = "~")
    )
    r2_d <- summary(stats::lm(reg2, data=d1))$r.squared
  }


  R2DXj.Xmj <- (r1_d - r2_d)/(1-r2_d)

  # --- f2DXj.Xmj
  f2DXj.Xmj <- R2DXj.Xmj/(1-R2DXj.Xmj)

  # --- R2YXj.XmjD
  reg1 <- stats::as.formula(
    paste(outcome, paste(c(bnch_reg,other_reg,treatment), collapse= "+"), sep = "~")
  )
  r1_y <- summary(stats::lm(reg1, data=d1))$r.squared

  if(is.null(other_reg)){
    r2_y <- 0
  } else{
    reg2 <- stats::as.formula(
      paste(outcome, paste(c(other_reg,treatment), collapse= "+"), sep = "~")
    )
    r2_y <- summary(stats::lm(reg2, data=d1))$r.squared

  }

  R2YXj.XmjD <- (r1_y - r2_y)/(1-r2_y)

  # --- f2YXj.XmjD
  f2YXj.XmjD <- R2YXj.XmjD/(1-R2YXj.XmjD)

  # --- fkD and kY [depends on kD]
  # Number of grid points to compute kYmax
  N <- 1000
  kD <- seq(from=0.1,to=((1-R2DXj.Xmj)/R2DXj.Xmj),length.out=N)
  kY <- numeric(N)
  for(i in 1:N){
    # Compute fkD
    fkD <- sqrt((kD[i]*R2DXj.Xmj)/(1-kD[i]*R2DXj.Xmj))
    # Compute kY
    kY[i] <- ((sqrt(1-(fkD^2)*f2DXj.Xmj))/sqrt(f2YXj.XmjD) - fkD*sqrt(f2DXj.Xmj))^2
  }


  # Return values
  return(
    data.frame(
      # Max kD
      kd_high=((1-R2DXj.Xmj)/R2DXj.Xmj),
      # Max kY
      ky_high=(max(kY,na.rm = TRUE))
    )
  )

}
