#' probability of conclusion being overturned using partial R2-based analysis with conditioning on treatment
#'
#' @param alpha significance level for hypothesis test (e.g. 0.05)
#' @param data data frame for analysis
#' @param outcome name of outcome variable
#' @param treatment name of treatment variable
#' @param bnch_reg name(s) of benchmark covariate(s)
#' @param other_reg name(s) of other covariate(s)
#' @param N number of points on grid = N^2
#' @param maxkd max of sensitivity parameter kD
#' @param maxky max of sensitivity parameter kY
#' @param k_kd mode (and median) of sensitivity parameter kD
#' @param k_ky mode (and median) of sensitivity parameter kY
#'
#' @import dplyr
#' @import tidyr
#'
#' @returns list with the following elements:
#' \item{dataplot}{data set used for contour plot}
#' \item{kdmax}{max of sensitivity parameter kD}
#' \item{kymax}{max of sensitivity parameter kY}
#' \item{frac_prob}{prob of conclusion being overturned (unwt)}
#' \item{frac_prob_wt}{prob of conclusion being overturned (wt)}
#' \item{frac_prob_rest}{prob of conclusion being overturned (unwt, rest)}
#' \item{frac_prob_rest_wt}{prob of conclusion being overturned (wt, rest)}
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
#' res4 <- salpr2cd(alpha=0.05,data=darfur,outcome=Y,treatment=D,bnch_reg=X,other_reg=X_oth,N=500)
#'
salpr2cd <- function(
    alpha,data,outcome,treatment,
    bnch_reg,other_reg,N,maxkd=NULL,
    maxky=NULL,k_kd=1,k_ky=1){


  # --- Compute ranges for kd,ky
  r1 <- ovbsa::kdkyrngpr2cd(
    data = data,
    outcome = outcome,
    treatment = treatment,
    bnch_reg = bnch_reg,
    other_reg = other_reg
  )

  skd_temp <- seq(from=0, to=r1$kd_high,length.out=N)
  skd <- skd_temp

  sky_temp <- seq(from=0, to=r1$ky_high,length.out=N)
  sky <- sky_temp

  # check that choice of k_kd is permissible
  if(k_kd<min(skd)|k_kd>max(skd)){
    stop("k_kd is outside permissible range")
  } else if(k_kd>(max(skd)/sqrt(exp(1)))){
    stop("k_kd is too large")
  }

  # check that choice of k_ky is permissible
  if(k_ky<min(sky)|k_ky>max(sky)){
    stop("k_ky is outside permissible range")
  } else if(k_ky>(max(sky)/sqrt(exp(1)))){
    stop("k_ky is too large")
  }

  # --------------- Estimate restricted model
  # Set up model
  mod1 <- stats::as.formula(
    paste(
      outcome, paste(c(treatment,bnch_reg, other_reg), collapse= "+"), sep = "~"
    )
  )
  # Results of regression
  res1 <- stats::lm(mod1,data = data)

  # Estimate
  est1 <- lmtest::coeftest(res1)[paste0(treatment),"Estimate"]


  # -------------- Generate data set for contour plots
  # To bind these variables locally set them initially as NULL
  X1 <- X2 <- Z1 <- Z2 <- Z3 <- Z4 <- NULL
  # Data set for contour plots
  d2 <- expand.grid(X1 = skd, X2 = sky) %>%
    mutate(
      # Adjusted estimate if unadjusted estimate>0
      Z1 = ovbsa::pr2cdbias(
        X1, X2, data = data, outcome = outcome,
        treatment = treatment, bnch_reg = bnch_reg,
        other_reg = other_reg, alpha = alpha
      )$adjestp,
      # Adjusted estimate if unadjusted estimate<0
      Z2 = ovbsa::pr2cdbias(
        X1, X2, data = data, outcome = outcome,
        treatment = treatment, bnch_reg = bnch_reg,
        other_reg = other_reg, alpha = alpha
      )$adjestn,
      # Lower boundary of CI if unadjusted estimate>0
      Z3 = ovbsa::pr2cdbias(
        X1, X2, data = data, outcome = outcome,
        treatment = treatment, bnch_reg = bnch_reg,
        other_reg = other_reg, alpha = alpha
      )$cilbp,
      # # Upper boundary of CI if unadjusted estimate<0
      Z4 = ovbsa::pr2cdbias(
        X1, X2, data = data, outcome = outcome,
        treatment = treatment, bnch_reg = bnch_reg,
        other_reg = other_reg, alpha = alpha
      )$ciubn
    ) %>%
    mutate(
      # est1>0: 1 if conf int lies to the right of zero, 0 otherwise
      cip_zero = ifelse((Z3>0),1,0),
      # est1<0: 1 if conf int lies to the left of zero, 0 otherwise
      cin_zero = ifelse((Z4<0),1,0)
    ) %>%
    mutate(
      linvwt_kd = ovbsa::linvx(x=X1, xvec=skd, k=k_kd),
      linvwt_ky = ovbsa::linvx(x=X2, xvec=sky, k=k_ky)
    ) %>%
    as.data.frame()



  # ----- Compute probability of conclusion of ---------------- #
  # ----- study being overturned: Unweighted  ----------------- #

  # ----- Case 1: No restriction on area of the contour plot ------- #

  # --- Unweighted
  # est1>0: fraction of area where confidence interval to the right of zero
  l1 <- d2[,"cip_zero"]
  l1_1 <- l1[!is.na(l1)]
  frac_pos <- sum(l1_1)/length(l1_1)

  # est1<0: fraction of area where confidence interval to the left of zero
  l2 <- d2[,"cin_zero"]
  l2_1 <- l2[!is.na(l2)]
  frac_neg <- sum(l2_1)/length(l2_1)

  fprob <- ifelse(est1>0,frac_pos,frac_neg)

  # ----- Weighted (Linear + Inverse)
  # est1>0: fraction of area where confidence interval to the right of zero
  l1 <- d2[,c("cip_zero","linvwt_kd","linvwt_ky")]
  l1$wt_kd <- l1[,c("linvwt_kd")]/sum(unique(l1[,c("linvwt_kd")]))
  l1$wt_ky <- l1[,c("linvwt_ky")]/sum(unique(l1[,c("linvwt_ky")]))
  l1$wt_kdky <- l1$wt_kd*l1$wt_ky
  l1_1 <- drop_na(l1)
  frac_pos_wt <- sum(l1_1[,c("cip_zero")]*l1_1[,c("wt_kdky")])/sum(l1_1[,c("wt_kdky")])

  # est1<0: fraction of area where confidence interval to the left of zero
  l2 <- d2[,c("cin_zero","linvwt_kd","linvwt_ky")]
  l2$wt_kd <- l2[,c("linvwt_kd")]/sum(unique(l1[,c("linvwt_kd")]))
  l2$wt_ky <- l2[,c("linvwt_ky")]/sum(unique(l1[,c("linvwt_ky")]))
  l2$wt_kdky <- l2$wt_kd*l2$wt_ky
  l2_1 <- drop_na(l2)
  frac_neg_wt <- sum(l2_1[,c("cin_zero")]*l2_1[,c("wt_kdky")])/sum(l2_1[,c("wt_kdky")])

  # weight
  fprob_wt <- ifelse(est1>0,frac_pos_wt,frac_neg_wt)



  # ----- Case 2: Restricted area of the contour plot ------- #

  if(!is.null(maxkd)&!is.null(maxky)){

    # check that choice of k_kd is permissible with restriction on maxkd
    if(k_kd<min(skd)|k_kd>maxkd){
      stop("k_kd is outside permissible range")
    } else if(k_kd>maxkd/sqrt(exp(1))){
      stop("k_kd is too large with restricted maxkD")
    }

    # check that choice of k_ky is permissible with restriction on maxky
    if(k_ky<min(skd)|k_ky>maxky){
      stop("k_ky is outside permissible range")
    } else if(k_ky>maxky/sqrt(exp(1))){
      stop("k_ky is too large with restricted maxkY")
    }

    # Construct restricted data frame
    d3 <- d2 %>%
      dplyr::filter(X1<=maxkd) %>%
      dplyr::filter(X2<=maxky) %>%
      as.data.frame()

    # ---- Unweighted
    # est1>0: fraction of area where confidence interval to the right of zero
    l1 <- d3[,"cip_zero"]
    l1_1 <- l1[!is.na(l1)]
    frac_pos <- sum(l1_1)/length(l1_1)

    # est1<0: fraction of area where confidence interval to the left of zero
    l2 <- d3[,"cin_zero"]
    l2_1 <- l2[!is.na(l2)]
    frac_neg <- sum(l2_1)/length(l2_1)

    fprob_rest <- ifelse(est1>0,frac_pos,frac_neg)

    # ----- Weighted (Linear + Inverse)
    # est1>0: fraction of area where confidence interval to the right of zero
    l1 <- d3[,c("cip_zero","linvwt_kd","linvwt_ky")]
    l1$wt_kd <- l1[,c("linvwt_kd")]/sum(unique(l1[,c("linvwt_kd")]))
    l1$wt_ky <- l1[,c("linvwt_ky")]/sum(unique(l1[,c("linvwt_ky")]))
    l1$wt_kdky <- l1$wt_kd*l1$wt_ky
    l1_1 <- drop_na(l1)
    frac_pos_wt <- sum(l1_1[,c("cip_zero")]*l1_1[,c("wt_kdky")])/sum(l1_1[,c("wt_kdky")])

    # est1<0: fraction of area where confidence interval to the left of zero
    l2 <- d3[,c("cin_zero","linvwt_kd","linvwt_ky")]
    l2$wt_kd <- l2[,c("linvwt_kd")]/sum(unique(l1[,c("linvwt_kd")]))
    l2$wt_ky <- l2[,c("linvwt_ky")]/sum(unique(l1[,c("linvwt_ky")]))
    l2$wt_kdky <- l2$wt_kd*l2$wt_ky
    l2_1 <- drop_na(l2)
    frac_neg_wt <- sum(l2_1[,c("cin_zero")]*l2_1[,c("wt_kdky")])/sum(l2_1[,c("wt_kdky")])

    fprob_rest_wt <- ifelse(est1>0,frac_pos_wt,frac_neg_wt)

  } else {
    fprob_rest <- NA
    fprob_rest_wt <- NA
  }




  # ---------- Return data and plots
  return(
    list(
      # data set used for plot
      dataplot = d2,
      # Max kD
      kdmax = r1$kd_high,
      # Max kY
      kymax = r1$ky_high,
      # Unweighted: How likely is OVB to make the true estimate zero?
      frac_prob = (1-fprob),
      # Weighted: How likely is OVB to make the true estimate zero?
      frac_prob_wt = (1-fprob_wt),
      # Unweighted: How likely is OVB to make the true estimate zero with restrictions on kD, kY?
      frac_prob_rest = (1-fprob_rest),
      # Weighted: How likely is OVB to make the true estimate zero with restrictions on kD, kY?
      frac_prob_rest_wt = (1-fprob_rest_wt)
    )
  )

}

