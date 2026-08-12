#' Compute bias adjusted confidence interval using truncated exponential prior distributions for kD, kY.
#'
#'
#' @param fit An object of class \code{lm}.
#' @param treatment A character string naming the treatment variable.
#' @param benchmark A character string naming the benchmark variable.
#' @param N Numeric value for grid size.
#' @param alpha Significance level.
#' @param medkd Median of the distribution of kD (default=1).
#' @param medky Median of the distribution of kY (default=1).
#'
#' @return A list containing the following two elements:
#' \item{results}{A data frame with 2 rows ("Unadjusted" and "Bias-adjusted") and
#' 2 columns ("Lower" and "Upper") containing the computed 100*(1-alpha)\%
#' unadjusted and bias-adjusted confidence intervals.}
#' \item{undstats}{A data frame of underlying statistics: the estimate, std error,
#' max(kD) and max(k(Y).}
#' \item{support_kdky_plot}{A \code{ggplot2} plot object visualizing the support of
#' the joint distribution of (kD,kY).}
#'
#'
#' @export
#'
#'@importFrom stats qnorm pnorm model.frame model.matrix model.response
#'
#' @examples
#' # CRAN requires checking if the optional package is installed before running the example
#' if (requireNamespace("sensemakr", quietly = TRUE)) {
#'   library(sensemakr)
#'
#'   # Define example variables
#'   myN <- 500
#'
#'   # Fit the model using the darfur dataset provided by sensemakr
#'   fit <- lm(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
#'               pastvoted + hhsize_darfur + female + village,
#'             data = sensemakr::darfur)
#'
#'   # Run baci function and see results
#'   res1 <- baci(fit = fit, treatment = "directlyharmed",
#'                benchmark = "female", N = myN, alpha = 1/100)
#'}
baci <- function(fit,
         treatment,
         benchmark,
         N = 1000, alpha = 5/100,
         medkd = 1, medky = 1){


  tictoc::tic("Extracting regression quantities ", "\n")
  # ---------- extract regression quantities ----------
  sm <- summary(fit)
  coefs <- sm$coefficients

  tau_hat <- coefs[treatment, "Estimate"]
  se_hat  <- coefs[treatment, "Std. Error"]
  df      <- sm$df[2]

  z <- qnorm(1 - alpha/2)

  # ---------- model frame / matrix ----------
  mf <- model.frame(fit)
  y  <- model.response(mf)
  X  <- model.matrix(fit)

  # treatment variable (D)
  d_col <- which(colnames(X) == treatment)

  # ---------- partial R2 of treatment with benchmark ----------
  # --------- r2_dx_j|xmj
  X_controls <- X[, -d_col, drop = FALSE]

  fit_d1 <- stats::lm(X[, d_col] ~ X_controls[, colnames(X_controls) != benchmark, drop = FALSE])
  fit_d2 <- stats::lm(X[, d_col] ~ X_controls)

  r2_d1 <- summary(fit_d1)$r.squared
  r2_d2 <- summary(fit_d2)$r.squared

  # using eqn 17 in CH (2020)
  r2d_bench <- (r2_d2 - r2_d1) / (1 - r2_d1)

  # ---------- partial R2 of outcome with benchmark ----------
  # --------- r2_yx_j|xmj
  cols_no_b <- dplyr::setdiff(colnames(X), benchmark)

  fit_y1 <- stats::lm(y ~ X[, cols_no_b, drop = FALSE])
  fit_y2 <- stats::lm(y ~ X)

  r2_y1 <- summary(fit_y1)$r.squared
  r2_y2 <- summary(fit_y2)$r.squared

  # using eqn 17 in CH (2020)
  r2y_bench <- (r2_y2 - r2_y1) / (1 - r2_y1)

  # ---------- partial R2 of outcome with treatment ----------
  # --------- r2_yd|x
  cols_no_d <- dplyr::setdiff(colnames(X), treatment)

  fit_t1 <- stats::lm(y ~ X[, cols_no_d, drop = FALSE])
  fit_t2 <- stats::lm(y ~ X)

  r2_t1 <- summary(fit_t1)$r.squared
  r2_t2 <- summary(fit_t2)$r.squared

  r2_yd_x <- (r2_t2 - r2_t1) / (1 - r2_t1)

  tictoc::toc()

  tictoc::tic("Computing values on the grid ", "\n")
  # ---------- grid ----------

  # kdmax: eqn 12 in Basu (2026)
  kdmax <- (1 - r2d_bench)/(r2d_bench)

  # kymax: eqn 12 in Basu (2026)
  kymax <- (1 - r2y_bench)/(r2y_bench)


  # Check if user supplied median(kD) and median(kY) are permissible
  if (!is.numeric(medkd) || medkd < 0 || medkd > kdmax) {
    stop(paste0("medkd must be a numeric value between 0 and ", kdmax, "."))
  }

  if (!is.numeric(medky) || medky < 0 || medky > kymax) {
    stop(paste0("medky must be a numeric value between 0 and ", kymax, "."))
  }

  # construct grid
  kd_grid <- seq(0, kdmax, length = N)
  ky_grid <- seq(0, kymax, length = N)
  grid <- expand.grid(kD = kd_grid, kY = ky_grid)

  # --- Here and after: u = z (in paper): unobserved confounder

  # compute r2_du|x: eqn 22 in CH (2020)
  grid$r2_du_x <- grid$kD * (r2d_bench/(1 - r2d_bench))

  # compute , r2_yu|x: eqn 71 in Online supplement, CH (2020)
  grid$r2_yu_x <- grid$kY * (r2y_bench/(1 - r2y_bench))

  # compute r2_yd|x (constant over kD,kY)
  grid$r2_yd_x <- r2_yd_x

  # compute r2_yu|dx: eqn 72 in Online supplement, CH (2020)
  grid$r2_yu_dx <-
    (sqrt(grid$r2_yu_x) - sqrt(grid$r2_yd_x) * sqrt(grid$r2_du_x))^2 /
    ((1 - grid$r2_yd_x) * (1 - grid$r2_du_x))

  # Define indicator for admissibility condition
  grid$is_admissible <- with(
    grid,
    r2_du_x > 0.01 & r2_du_x < 0.99 &
      r2_yu_x > 0.01 & r2_yu_x < 0.99 &
      r2_yu_dx > 0.01 & r2_yu_dx < 0.99
  )

  # ---------- Create plot for admissible grid boundary ----------
  plot_grid <- grid
  plot_grid$admissible_num <- as.numeric(plot_grid$is_admissible)

  # Centroid calculation to place interior text dynamically
  admissible_pts <- subset(plot_grid, is_admissible)
  mean_kD <- mean(admissible_pts$kD)
  mean_kY <- mean(admissible_pts$kY)

  support_kdky_plot <- ggplot2::ggplot(plot_grid, ggplot2::aes(x = kD, y = kY)) +
    ggplot2::geom_tile(ggplot2::aes(fill = is_admissible), alpha = 0.25) +
    ggplot2::geom_contour(
      ggplot2::aes(z = admissible_num),
      breaks = 0.5,
      color = "navy blue",
      linewidth = 1
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "skyblue", "FALSE" = "white"),
      guide = "none"
    ) +
    ggplot2::annotate(
      "text",
      x = mean_kD,
      y = mean_kY,
      label = "Admissible Grid Region",
      fontface = "bold",
      color = "navy blue",
      size = 4.5
    ) +
    ggplot2::labs(
      title = "Admissible Grid Boundary in (k_D, k_Y) Space",
      x = expression(k[D]),
      y = expression(k[Y])
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # retain admissible grid points for CI calculation
  grid <- subset(grid, is_admissible)

  # ---------- apply bias formula on filtered grid ----------
  # absolute value of bias: eqn 13, CH (2020)
  bias_mag <- se_hat * sqrt(df) *
    sqrt((grid$r2_yu_dx * grid$r2_du_x) / (1 - grid$r2_du_x))

  # choose worst-case sign
  bias <- sign(tau_hat) * bias_mag

  # compute bias adjusted estimate
  est_adj <- tau_hat - bias

  # compute bias adjusted std error: eqn 12, CH (2020)
  se_adj <- se_hat * sqrt(df/(df - 1)) *
    sqrt((1 - grid$r2_yu_dx) / (1 - grid$r2_du_x))

  # compute lower boundary of CI: eqn 10, 11 in Basu (2026)
  lower <- est_adj - z * se_adj

  # compute upper boundary of CI: eqn 10, 11 in Basu (2026)
  upper <- est_adj + z * se_adj

  # initial output
  out <- data.frame(
    kD = grid$kD,
    kY = grid$kY,
    r2_du_x = grid$r2_du_x,
    r2_yu_x = grid$r2_yu_x,
    r2_yd_x = grid$r2_yd_x,
    r2_yu_dx = grid$r2_yu_dx,
    estimate = est_adj,
    se = se_adj,
    lower = lower,
    upper = upper
  )
  tictoc::toc()


  # --- prior-weighted mixture CDFs over the (kD, kY) grid
  # --- with different prior distributions for kD and kY



  # --- truncated exponendtial distributions: exp(1)
  tictoc::tic("Prior distribution: truncated exponential ", "\n")
  # --- see section 3.2, Basu (2026)
  # solve for beta in the distribution for kD
  mybetaexp_kd <- stats::uniroot(
    function(x) 1 - exp(-(medkd/x)) - pexp(kdmax)/2,
    lower = 0.01,
    upper= kdmax
  )$root

  # solve for beta in the distribution for kY
  mybetaexp_ky <- stats::uniroot(
    function(x) 1 - exp(-(medky/x)) - pexp(kymax)/2,
    lower = 0.01,
    upper= kymax
  )$root

  # weights
  out$w_e <- dtruncexp(out$kD, rate=mybetaexp_kd) * dtruncexp(out$kY, rate=mybetaexp_ky)
  out$w_e <- out$w_e / sum(out$w_e)

  # Root finding
  prior_lower_e <- stats::uniroot(
    # eqn 19, Basu (2026)
    function(x) sum(out$w_e * pnorm((x - out$estimate)/out$se)) - (alpha/2),
    lower = (min(out$lower)-5),
    upper= (max(out$upper)+5)
  )$root

  prior_upper_e <- stats::uniroot(
    # eqn 19, Basu (2026)
    function(x) sum(out$w_e * pnorm((x - out$estimate)/out$se)) - (1-(alpha/2)),
    lower = (min(out$lower)-5),
    upper= (max(out$upper)+5)
  )$root

  tictoc::toc()



  # -- 100*(1-alpha) bias adjusted confidence interval
  ci_adj <- c(prior_lower_e, prior_upper_e)

  # -- 100*(1-alpha) bias adjusted confidence interval
  ci_unadj <- c(
    tau_hat - abs(stats::qt((alpha/2),df=df,lower.tail = TRUE))*se_hat,
    tau_hat + abs(stats::qt((alpha/2),df=df,lower.tail = TRUE))*se_hat
    )


  biasadjci <- as.data.frame(
    rbind(ci_unadj,ci_adj)
  )
  colnames(biasadjci) <- c("Lower","Upper")
  rownames(biasadjci) <- c("Unadjusted CI","Bias-adjusted CI")

  # --- Underlying stats
  und_stats <- as.data.frame(c(tau_hat,se_hat,kdmax,kymax))
  rownames(und_stats) <- c(
    "estimate","std error","max(kD)","max(kY)"
    )
  colnames(und_stats) <- c("Values")


  cat("Results available now!", "\n")

  # Return named list
  return(list(
    results = biasadjci,
    undstats = und_stats,
    support_kdky_plot = support_kdky_plot
  ))

}
