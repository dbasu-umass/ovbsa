#' bias and std error for (kd,ky) using partial R2-based analysis with conditioning on treatment
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
#' require("sensemakr")
#' Y <- "peacefactor"
#' D <- "directlyharmed"
#' X <- "female"
#' X_oth <- c("village","age","farmer_dar","herder_dar","pastvoted","hhsize_darfur")
#'
#'
#' res5<-pr2cdbias(kd=1,ky=1,alpha=0.05,data=darfur,outcome=Y,treatment=D,bnch_reg=X,other_reg=X_oth)
#'
pr2cdbias <- function(
    kd,ky,alpha,data,outcome,treatment,
    bnch_reg,other_reg=NULL){


  # ----------------------
  # ------ Create data set
  d1 <- as.data.frame(data)

  # -----------------------------------------
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


  # ----- Case 3: partial r2-based benchmarking
  # ----- conditioning on D (treatment)

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

  # Use definition [equation 17]
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

  # Use definition [equation 17]
  R2YXj.Xmj <- (r1_y - r2_y)/(1-r2_y)



  # ---- R2YD_X

  # First residual: Y_X
  YX <- stats::as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  u1 <- stats::lm(YX, data = d1)$residuals
  # remove NAs
  u1 <- u1[!is.na(u1)]

  # Second residual: D_X
  DX <- stats::as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  u2 <- stats::lm(DX, data = d1)$residuals
  # remove NAs
  u2 <- u2[!is.na(u2)]

  # compute R2YD_X
  N <- min(length(u1),length(u2))
  R2YD_X <- summary(stats::lm(u1[1:N]~u2[1:N]))$r.squared


  # ---- Compute R2DZ_X
  R2DZ_X <- kd*(R2DXj.Xmj/(1-R2DXj.Xmj))

  # ---- Compute R2YZ_X
  R2YZ_X <- ky*(R2YXj.Xmj/(1-R2YXj.Xmj))

  # ----- Compute R2YZ_DX
  # fkd
  fkd <- (sqrt(kd*R2DXj.Xmj))/(sqrt(1-kd*R2DXj.Xmj))

  # R2YXj.[Xmj,D]
  reg1 <- stats::as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg, treatment), collapse= "+"), sep = "~")
  )
  r1_y <- summary(stats::lm(reg1, data=d1))$r.squared

  reg2 <- stats::as.formula(
    paste(outcome, paste(c(other_reg,treatment), collapse= "+"), sep = "~")
  )
  r2_y <- summary(stats::lm(reg2, data=d1))$r.squared

  R2YXj.XmjD <- (r1_y - r2_y)/(1-r2_y)

  # fDXj.Xmj
  f2DXj.Xmj <- sqrt(R2DXj.Xmj/(1-R2DXj.Xmj))

  # eta
  myeta <- (sqrt(ky)+abs(fkd*f2DXj.Xmj))/(sqrt(1-(fkd*f2DXj.Xmj)^2))

  # R2YZ_DX
  R2YZ_DX <- (myeta^2)*(R2YXj.XmjD/(1-R2YXj.XmjD))


  # --- To be used for computing bias and adjusted estimate
  R2DZ_X_3 <- R2DZ_X
  R2YZ_DX_3 <- R2YZ_DX


  # -----------------------------------------
  # ---------- Compute bias and se of adj est

  # --- Absolute value of omitted variable bias
  bias_abs_3 <- (se1)*sqrt((dof1*R2YZ_DX_3*R2DZ_X_3)/(1-R2DZ_X_3))

  # tau_hat if est1>0
  tau_hat_3_pos <- (est1 - bias_abs_3)
  # tau_hat if est1<0
  tau_hat_3_neg <- (est1 + bias_abs_3)

  # Std Err of tau_hat
  A <- (1-R2YZ_DX_3)/(1-R2DZ_X_3)
  B <- dof1/(dof1-1)
  suppressWarnings(se_tauhat_3 <- se1*sqrt(A*B))

  # Lower bound of CI if est1>0
  ci_lb_p <- tau_hat_3_pos - abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_3

  # Upper bound of CI if est1>0
  ci_ub_p <- tau_hat_3_pos + abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_3

  # Lower bound of CI if est1<0
  ci_lb_n <- tau_hat_3_neg - abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_3

  # Upper bound of CI if est1<0
  ci_ub_n <- tau_hat_3_neg + abs(stats::qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_3

  # Return
  return(
    list(
      # Adj std err when unad est>0
      adjestp=tau_hat_3_pos,
      # Adj std err when unad est<0
      adjestn=tau_hat_3_neg,
      # Adj lower boundary of CI when unad est>0
      cilbp=ci_lb_p,
      # Adj upper boundary of CI when unad est>0
      ciubp=ci_ub_p,
      # Adj lower boundary of CI when unad est<0
      cilbn=ci_lb_n,
      # Adj upper boundary of CI when unad est<0
      ciubn=ci_ub_n
    )
  )

}

