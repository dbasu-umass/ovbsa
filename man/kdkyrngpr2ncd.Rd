% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kdkyrngpr2ncd.R
\name{kdkyrngpr2ncd}
\alias{kdkyrngpr2ncd}
\title{compute max(kD) and max(kY) for partial R2-based analysis without conditioning on treatment}
\usage{
kdkyrngpr2ncd(data, outcome, treatment, bnch_reg, other_reg = NULL)
}
\arguments{
\item{data}{data frame for analysis}

\item{outcome}{name of outcome variable}

\item{treatment}{name of treatment variable}

\item{bnch_reg}{name(s) of benchmark covariate(s)}

\item{other_reg}{name(s) of other covariates}
}
\value{
a data frame with 2 columns and 1 row:
\item{kd_high}{max(kD), a scalar}
\item{ky_high}{max(kY), a scalar}
}
\description{
compute max(kD) and max(kY) for partial R2-based analysis without conditioning on treatment
}
\examples{
require("sensemakr")
Y <- "peacefactor"
D <- "directlyharmed"
X <- "female"
X_oth <- c("village","age","farmer_dar","herder_dar","pastvoted","hhsize_darfur")


r1 <- kdkyrngpr2ncd(data=darfur,outcome=Y,treatment=D,bnch_reg=X,other_reg=X_oth)

}
