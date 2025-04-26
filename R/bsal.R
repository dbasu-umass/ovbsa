bsal <- function(kd,ky,alpha,
                   data,outcome,treatment,
                   bnch_reg, other_reg){


  # ------ Create data set
  d1 <- as.data.frame(data)

  # ------ Estimate the restricted regression
  # Set up model
  mod1 <- as.formula(
    paste(
      outcome, paste(c(treatment,bnch_reg, other_reg), collapse= "+"), sep = "~"
    )
  )
  # Results of regression
  res1 <- lm(mod1,data = d1)

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
  DXj <- as.formula(paste(treatment, paste(bnch_reg, collapse= "+"), sep = "~"))
  R2DXj <- summary(lm(DXj, data = d1))$r.squared

  # ---- R2DX
  DX <- as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  R2DX <- summary(lm(DX,data = d1))$r.squared

  # ---- R2YXj
  YXj <- as.formula(paste(outcome, paste(bnch_reg, collapse= "+"), sep = "~"))
  R2YXj <- summary(lm(YXj, data = d1))$r.squared


  # ----- R2YX
  YX <- as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  R2YX <- summary(lm(YX,data = d1))$r.squared

  # ---- R2YD_X
  # First residual: Y_X
  u1 <- lm(YX, data = d1)$residuals
  u1 <- u1[!is.na(u1)]
  # Second residual: D_X
  u2 <- lm(DX, data = d1)$residuals
  u2 <- u2[!is.na(u2)]
  # R2YD_X
  N <- min(length(u1),length(u2))
  R2YD_X <- summary(lm(u1[1:N]~u2[1:N]))$r.squared


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
  reg1 <- as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  r1_d <- summary(lm(reg1, data=d1))$r.squared

  if(is.null(other_reg)){
    r2_d <- 0
  } else{
    reg2 <- as.formula(
      paste(treatment, paste(other_reg, collapse= "+"), sep = "~")
    )
    r2_d <- summary(lm(reg2, data=d1))$r.squared
  }


  R2DXj.Xmj <- (r1_d - r2_d)/(1-r2_d)


  # ---- R2YXj.Xmj
  reg1 <- as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  r1_y <- summary(lm(reg1, data=d1))$r.squared

  if(is.null(other_reg)){
    r2_y <- 0
  } else{
    reg2 <- as.formula(
      paste(outcome, paste(c(other_reg), collapse= "+"), sep = "~")
    )
    r2_y <- summary(lm(reg2, data=d1))$r.squared

  }

  R2YXj.Xmj <- (r1_y - r2_y)/(1-r2_y)



  # ---- R2YD_X
  # First residual: Y_X
  YX <- as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  u1 <- lm(YX, data = d1)$residuals
  u1 <- u1[!is.na(u1)]
  # Second residual: D_X
  DX <- as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  u2 <- lm(DX, data = d1)$residuals
  u2 <- u2[!is.na(u2)]
  # R2YD_X
  N <- min(length(u1),length(u2))
  R2YD_X <- summary(lm(u1[1:N]~u2[1:N]))$r.squared


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



  # ----- Case 3: partial r2-based benchmarking
  # ----- conditioning on D (treatment)

  # --- R2DXj.Xmj
  reg1 <- as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  r1_d <- summary(lm(reg1, data=d1))$r.squared

  if(is.null(other_reg)){
    r2_d <- 0
  } else{
    reg2 <- as.formula(
      paste(treatment, paste(other_reg, collapse= "+"), sep = "~")
    )
    r2_d <- summary(lm(reg2, data=d1))$r.squared
  }

  # Use definition [equation 17]
  R2DXj.Xmj <- (r1_d - r2_d)/(1-r2_d)


  # ---- R2YXj.Xmj
  reg1 <- as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  r1_y <- summary(lm(reg1, data=d1))$r.squared

  if(is.null(other_reg)){
    r2_y <- 0
  } else{
    reg2 <- as.formula(
      paste(outcome, paste(c(other_reg), collapse= "+"), sep = "~")
    )
    r2_y <- summary(lm(reg2, data=d1))$r.squared

  }

  # Use definition [equation 17]
  R2YXj.Xmj <- (r1_y - r2_y)/(1-r2_y)



  # ---- R2YD_X

  # First residual: Y_X
  YX <- as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  u1 <- lm(YX, data = d1)$residuals
  # remove NAs
  u1 <- u1[!is.na(u1)]

  # Second residual: D_X
  DX <- as.formula(
    paste(treatment, paste(c(bnch_reg, other_reg), collapse= "+"), sep = "~")
  )
  u2 <- lm(DX, data = d1)$residuals
  # remove NAs
  u2 <- u2[!is.na(u2)]

  # compute R2YD_X
  N <- min(length(u1),length(u2))
  R2YD_X <- summary(lm(u1[1:N]~u2[1:N]))$r.squared


  # ---- Compute R2DZ_X
  R2DZ_X <- kd*(R2DXj.Xmj/(1-R2DXj.Xmj))
  if(R2DZ_X>1){
    stop("Case 3: R2DZ_X>1. Choose a lower kd.")
  }

  # ---- Compute R2YZ_X
  R2YZ_X <- ky*(R2YXj.Xmj/(1-R2YXj.Xmj))
  if(R2YZ_X>1){
    stop("Case 3: R2YZ_X>1. Choose a lower kd.")
  }

  # ----- Compute R2YZ_DX
  # fkd
  fkd <- (sqrt(kd*R2DXj.Xmj))/(sqrt(1-kd*R2DXj.Xmj))

  # R2YXj.[Xmj,D]
  reg1 <- as.formula(
    paste(outcome, paste(c(bnch_reg, other_reg, treatment), collapse= "+"), sep = "~")
  )
  r1_y <- summary(lm(reg1, data=d1))$r.squared

  reg2 <- as.formula(
    paste(outcome, paste(c(other_reg,treatment), collapse= "+"), sep = "~")
  )
  r2_y <- summary(lm(reg2, data=d1))$r.squared

  R2YXj.XmjD <- (r1_y - r2_y)/(1-r2_y)

  # fDXj.Xmj
  f2DXj.Xmj <- sqrt(R2DXj.Xmj/(1-R2DXj.Xmj))

  # eta
  myeta <- (sqrt(ky)+(fkd*f2DXj.Xmj))/(sqrt(1-(fkd*f2DXj.Xmj)^2))

  # R2YZ_DX
  R2YZ_DX <- (myeta^2)*(R2YXj.XmjD/(1-R2YXj.XmjD))
  if(R2YZ_DX>1){
    stop("Case 3: R2YZ_DX>1. Choose a lower kd and/or ky.")
  }


  # --- To be used for computing bias and adjusted estimate
  R2DZ_X_3 <- R2DZ_X
  R2YZ_DX_3 <- R2YZ_DX


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

  bias_abs_3 <- (se1)*sqrt((dof1*R2YZ_DX_3*R2DZ_X_3)/(1-R2DZ_X_3))
  bias_3 <- ifelse(
    est1>=0, bias_abs_3, (-1)*bias_abs_3
  )

  # ------ adjusted estimate = tau_hat
  tau_hat_1 <- (est1 - bias_1)
  tau_hat_2 <- (est1 - bias_2)
  tau_hat_3 <- (est1 - bias_3)

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

  if(R2YZ_DX_3>1){
    stop("Case 3: R2YZ|DX>1. Reduce kd and/or ky.")
  } else {
    A <- (1-R2YZ_DX_3)/(1-R2DZ_X_3)
    B <- dof1/(dof1-1)
    se_tauhat_3 <- se1*sqrt(A*B)
  }


  # ------ Bounds of conf int
  ci_lb_1 <- tau_hat_1 - abs(qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_1
  ci_ub_1 <- tau_hat_1 + abs(qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_1

  ci_lb_2 <- tau_hat_2 - abs(qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_2
  ci_ub_2 <- tau_hat_2 + abs(qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_2

  ci_lb_3 <- tau_hat_3 - abs(qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_3
  ci_ub_3 <- tau_hat_3 + abs(qt((alpha/2),df=dof1,lower.tail = TRUE))*se_tauhat_3


  # -------------------------------- #
  # ------------ Results ----------- #

  # Combined results
  comb_results <- cbind(
    c(R2YD_X,R2YD_X,R2YD_X),
    c(R2DZ_X_1,R2DZ_X_2,R2DZ_X_3),
    c(R2YZ_DX_1,R2YZ_DX_2,R2YZ_DX_3),
    c(est1,est1,est1),
    c(tau_hat_1,tau_hat_2,tau_hat_3),
    c(se_tauhat_1,se_tauhat_2,se_tauhat_3),
    c(ci_lb_1,ci_lb_2,ci_lb_3),
    c(ci_ub_1,ci_ub_2,ci_ub_3)
  )

  colnames(comb_results) <- c("r2yd.x","r2dz.x","r2yz.dx",
                              "estimate",
                              "adjusted_estimate","adjusted_se",
                              "adjusted_lower_CI",
                              "adjusted_upper_CI")
  rownames(comb_results) <- c("Case 1", "Case 2", "Case 3")

  # Return result as matrix
  return(t(comb_results))

}
