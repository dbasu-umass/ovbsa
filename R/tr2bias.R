#' bias and std error for (kd,ky) using total R2-based analysis
#'
#' @param kd sensitivity parameter kD (scalar)
#' @param ky sensitivity parameter kY (scalar)
#' @param alpha significance level for hypothesis test (e.g. 0.05)
#' @param data data frame for analysis
#' @param outcome name of outcome variable
#' @param treatment name of treatment variable
#' @param bnch_reg name(s) of benchmark covariate(s)
#' @param other_reg name(s) of other covariate(s)
#'
#' @returns a list with the following elements:
#' \item{adjestp}{Adj std error when unadj estimate>0}
#' \item{adjestn}{Adj std error when unadj estimate<0}
#' \item{cilbp}{Adj lower boundary of conf int when unadj estimate>0}
#' \item{ciubp}{Adj upper boundary of conf int when unadj estimate>0}
#' \item{cilbn}{Adj lower boundary of conf int when unadj estimate<0}
#' \item{ciubn}{Adj upper boundary of conf int when unadj estimate<0}
#' @export
#'
#' @examples
#'
#' if(require("sensemakr")){
#' Y <- "peacefactor"
#' D <- "directlyharmed"
#' X <- "female"
#' X_oth <- c("village","age","farmer_dar","herder_dar","pastvoted","hhsize_darfur")
#' }
#'
#' res2 <- tr2bias(kd=1,ky=1,alpha=0.05,data=darfur,outcome=Y,treatment=D,bnch_reg=X,other_reg=X_oth)
#'
tr2bias <- function(kd,ky,alpha,data,outcome,treatment,bnch_reg,other_reg=NULL){

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

  # ----- R2YZ_X
  R2YZ_X <- ky*(R2YXj/(1-R2YX))

  # ----- R2YZ_DX
  num1 <- (sqrt(R2YZ_X) - sqrt(R2YD_X)*sqrt(R2DZ_X))^2
  den1 <- (1-R2YD_X)*(1-R2DZ_X)
  R2YZ_DX <- num1/den1

  # --- To be used for computing bias and adjusted estimate
  R2DZ_X_1 <- R2DZ_X
  R2YZ_DX_1 <- R2YZ_DX



  # -----------------------------------------
  # ---------- Compute bias and se of adj est

  # Absolute bias
  bias_abs_1 <- (se1)*sqrt((dof1*R2YZ_DX_1*R2DZ_X_1)/(1-R2DZ_X_1))

  # tau_hat if est1>0
  tau_hat_1_pos <- (est1 - bias_abs_1)
  # tau_hat if est1<0
  tau_hat_1_neg <- (est1 + bias_abs_1)

  # Std Err of tau_hat
  A <- (1-R2YZ_DX_1)/(1-R2DZ_X_1)
  B <- dof1/(dof1-1)
  suppressWarnings(se_tauhat_1 <- se1*sqrt(A*B))

  # Lower bound of CI if est1>0
  ci_lb_p <- tau_hat_1_pos - abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_1

  # Upper bound of CI if est1>0
  ci_ub_p <- tau_hat_1_pos + abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_1

  # Lower bound of CI if est1<0
  ci_lb_n <- tau_hat_1_neg - abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_1

  # Upper bound of CI if est1<0
  ci_ub_n <- tau_hat_1_neg + abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_1

  # Return
  return(
    list(
      # Adj std err when unadj est>0
      adjestp=tau_hat_1_pos,
      # Adj std err when unadj est<0
      adjestn=tau_hat_1_neg,
      # Adj lower boundary of CI when unadj est>0
      cilbp=ci_lb_p,
      # Adj upper boundary of CI when unadj est>0
      ciubp=ci_ub_p,
      # Adj lower boundary of CI when unadj est<0
      cilbn=ci_lb_n,
      # Adj upper boundary of CI when unadj est<0
      ciubn=ci_ub_n
    )
  )

}

