#' basic sensitivity analysis of omitted variable bias
#'
#' @param kd sensitivity parameter kD (scalar)
#' @param ky sensitivity parameter kY (scalar)
#' @param alpha significance level for hypothesis test (e.g. 0.05)
#' @param data data frame for analysis
#' @param outcome name of outcome variable
#' @param treatment name of treatment variable
#' @param bnch_reg name(s) of benchmark covariate(s)
#' @param other_reg name(s) of other regressors
#'
#' @returns a matrix with following rows for case 1, 2 and 3 (in columns):
#' \item{r2yd.x}{partial R2 of Y on D conditioning on X}
#' \item{r2dz.x}{partial R2 of D on Z conditioning on X}
#' \item{r2yz.dx}{partial R2 of Y on Z conditioning on D and X}
#' \item{estimate}{unadjusted parameter estimate}
#' \item{adjusted_estimate}{bias-adjusted parameter estimate}
#' \item{adjusted_se}{bias-adjusted standard error}
#' \item{adjusted_lower_CI}{bias-adjusted confidence interval lower boundary}
#' \item{adjusted_upper_CI}{bias-adjusted confidence interval upper boundary}
#'
#' @export
#'
#' @examples
#'
#' require("sensemakr")
#' Y <- "peacefactor"
#' D <- "directlyharmed"
#' X <- "female"
#' X_oth <- c("village","age","farmer_dar","herder_dar","pastvoted","hhsize_darfur")
#'
#'
#' res1 <- bsal(kd=1,ky=1,alpha=0.05,data=darfur,outcome=Y,treatment=D,bnch_reg=X,other_reg=X_oth)
#'
bsal <- function(kd,ky,alpha,
                   data,outcome,treatment,
                   bnch_reg, other_reg){


  # ------ Create data set
  d1 <- as.data.frame(data)

  # ------ Estimate the restricted regression
  # Set up model
  mod1 <- stats::as.formula(
    paste(
      outcome, paste(c(treatment,bnch_reg, other_reg), collapse= "+"), sep = "~"
    )
  )
  # Results of regression
  res1 <- stats::lm(mod1,data = d1)

  # Estimate
  est1 <- lmtest::coeftest(res1)[paste0(treatment),"Estimate"]
  # Std Error
  se1 <- lmtest::coeftest(res1)[paste0(treatment),"Std. Error"]
  # T-statistic
  tstat1 <- lmtest::coeftest(res1)[paste0(treatment),"t value"]
  # Degrees of freedom
  dof1 <- res1$df.residual
  # Number of regressors
  num_reg <- length(res1$coefficients)

  # ---------------------------------------------- #
  # ---------- Partial r2 measures --------------- #


  # -------- Case 1: total r2-based benchmarking

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

  # ---- R2YD_X
  # First residual: Y_X
  u1 <- stats::lm(YX, data = d1)$residuals
  u1 <- u1[!is.na(u1)]
  # Second residual: D_X
  u2 <- stats::lm(DX, data = d1)$residuals
  u2 <- u2[!is.na(u2)]
  # R2YD_X
  N <- min(length(u1),length(u2))
  R2YD_X <- summary(stats::lm(u1[1:N]~u2[1:N]))$r.squared


  # ---- R2DZ_X
  R2DZ_X <- kd*(R2DXj/(1-R2DX))
  if(R2DZ_X>1){
    stop("Case 1: R2DZ_X>1. Choose a lower kd")
  }

  # ----- R2YZ_X
  R2YZ_X <- ky*(R2YXj/(1-R2YX))
  if(R2YZ_X>1){
    stop("Case 1: R2YZ_X>1. Choose a lower ky")
  }

  # ----- R2YZ_DX
  num1 <- (sqrt(R2YZ_X) - sqrt(R2YD_X)*sqrt(R2DZ_X))^2
  den1 <- (1-R2YD_X)*(1-R2DZ_X)
  R2YZ_DX <- num1/den1
  if(R2YZ_DX>1){
    stop("Case 1: R2YZ|DX>1. Choose a lower kd and/or ky.")
  }

  # --- To be used for computing bias and adjusted estimate
  R2DZ_X_1 <- R2DZ_X
  R2YZ_DX_1 <- R2YZ_DX



  # -------- Case 2: partial r2-based benchmarking
  # -------- not conditioning on D (treatment)

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


  # ---- R2YXj.Xmj
  reg1 <- stats::as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  r1_y <- summary(stats::lm(reg1, data=d1))$r.squared

  if(is.null(other_reg)){
    r2_y <- 0
  } else{
    reg2 <- stats::as.formula(
      paste(outcome, paste(c(other_reg), collapse= "+"), sep = "~")
    )
    r2_y <- summary(stats::lm(reg2, data=d1))$r.squared

  }

  R2YXj.Xmj <- (r1_y - r2_y)/(1-r2_y)



  # ---- R2YD_X
  # First residual: Y_X
  YX <- stats::as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  u1 <- stats::lm(YX, data = d1)$residuals
  u1 <- u1[!is.na(u1)]
  # Second residual: D_X
  DX <- stats::as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  u2 <- stats::lm(DX, data = d1)$residuals
  u2 <- u2[!is.na(u2)]
  # R2YD_X
  N <- min(length(u1),length(u2))
  R2YD_X <- summary(stats::lm(u1[1:N]~u2[1:N]))$r.squared


  # Compute R2DZ_X
  R2DZ_X <- kd*(R2DXj.Xmj/(1-R2DXj.Xmj))
  if(R2DZ_X>1){
    stop("Case 2: R2DZ_X>1. Choose a lower kd.")
  }

  # Compute R2YZ_X
  R2YZ_X <- ky*(R2YXj.Xmj/(1-R2YXj.Xmj))
  if(R2YZ_X>1){
    stop("Case 2: R2YZ_X>1. Choose a lower kd.")
  }

  # Compute R2YZ_DX
  num1 <- (sqrt(R2YZ_X) - sqrt(R2YD_X)*sqrt(R2DZ_X))^2
  den1 <- (1-R2YD_X)*(1-R2DZ_X)
  R2YZ_DX <- num1/den1
  if(R2YZ_DX>1){
    stop("Case 2: R2YZ|DX>1. Choose a lower kd and/or ky.")
  }

  # --- To be used for computing bias and adjusted estimate
  R2DZ_X_2 <- R2DZ_X
  R2YZ_DX_2 <- R2YZ_DX



  # ------------------------------------------------------- #
  # -------- Adjusted estimate, SE and bounds of CI ------- #

  # ---- Bias [equation 13]
  bias_abs_1 <- (se1)*sqrt((dof1*R2YZ_DX_1*R2DZ_X_1)/(1-R2DZ_X_1))
  bias_1 <- ifelse(
    est1>=0, bias_abs_1, (-1)*bias_abs_1
  )

  bias_abs_2 <- (se1)*sqrt((dof1*R2YZ_DX_2*R2DZ_X_2)/(1-R2DZ_X_2))
  bias_2 <- ifelse(
    est1>=0, bias_abs_2, (-1)*bias_abs_2
  )

  # ------ adjusted estimate = tau_hat
  tau_hat_1 <- (est1 - bias_1)
  tau_hat_2 <- (est1 - bias_2)
  #tau_hat_3 <- (est1 - bias_3)

  # ------ adjusted SE = se(tau_hat)
  if(R2YZ_DX_1>1){
    stop("Case 1: R2YZ|DX>1. Reduce kd and/or ky.")
  } else {
    A <- (1-R2YZ_DX_1)/(1-R2DZ_X_1)
    B <- dof1/(dof1-1)
    se_tauhat_1 <- se1*sqrt(A*B)
  }

  if(R2YZ_DX_2>1){
    stop("Case 2: R2YZ|DX>1. Reduce kd and/or ky.")
  } else {
    A <- (1-R2YZ_DX_2)/(1-R2DZ_X_2)
    B <- dof1/(dof1-1)
    se_tauhat_2 <- se1*sqrt(A*B)
  }


  # ------ Bounds of conf int
  ci_lb_1 <- tau_hat_1 - abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_1
  ci_ub_1 <- tau_hat_1 + abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_1

  ci_lb_2 <- tau_hat_2 - abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_2
  ci_ub_2 <- tau_hat_2 + abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_2


  # -------------------------------- #
  # ------------ Results ----------- #


  # Combined results
  comb_results <- cbind(
    # row 1
    c(R2YD_X,R2YD_X),
    # row 2
    c(R2DZ_X_1,R2DZ_X_2),
    # row 3
    c(R2YZ_DX_1,R2YZ_DX_2),
    # row 4
    c(est1,est1),
    # row 5
    c(tau_hat_1,tau_hat_2),
    # row 6
    c(se_tauhat_1,se_tauhat_2),
    # row 7
    c(ci_lb_1,ci_lb_2),
    # row 8
    c(ci_ub_1,ci_ub_2)
  )

  colnames(comb_results) <- c("r2yd.x","r2dz.x","r2yz.dx",
                              "estimate",
                              "adjusted_estimate","adjusted_se",
                              "adjusted_lower_CI",
                              "adjusted_upper_CI")
  rownames(comb_results) <- c("Case 1", "Case 2")

  # Return result as matrix
  return(t(comb_results))

}
