% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saltr2.R
\name{saltr2}
\alias{saltr2}
\title{probability of conclusion being overturned using total R2-based analysis}
\usage{
saltr2(
  alpha,
  data,
  outcome,
  treatment,
  bnch_reg,
  other_reg,
  N,
  maxkd = NULL,
  maxky = NULL,
  k_kd = 1,
  k_ky = 1
)
}
\arguments{
\item{alpha}{significance level for hypothesis test (e.g. 0.05)}

\item{data}{data frame for analysis}

\item{outcome}{name of outcome variable}

\item{treatment}{name of treatment variable}

\item{bnch_reg}{name(s) of benchmark covariate(s)}

\item{other_reg}{name(s) of other covariate(s)}

\item{N}{number of points on grid = N^2}

\item{maxkd}{max of sensitivity parameter kD}

\item{maxky}{max of sensitivity parameter kY}

\item{k_kd}{mode (and median) of sensitivity parameter kD}

\item{k_ky}{mode (and median) of sensitivity parameter kY}
}
\value{
list with the following elements:
\item{dataplot}{data set used for contour plot}
\item{kdmax}{max of sensitivity parameter kD}
\item{kymax}{max of sensitivity parameter kY}
\item{frac_prob}{prob of conclusion being overturned (unwt)}
\item{frac_prob_wt}{prob of conclusion being overturned (wt)}
\item{frac_prob_rest}{prob of conclusion being overturned (unwt, rest)}
\item{frac_prob_rest_wt}{prob of conclusion being overturned (wt, rest)}
}
\description{
probability of conclusion being overturned using total R2-based analysis
}
\examples{

require("sensemakr")
Y <- "peacefactor"
D <- "directlyharmed"
X <- "female"
X_oth <- c("village","age","farmer_dar","herder_dar","pastvoted","hhsize_darfur")

darfur1 <- dplyr::slice_sample(darfur, prop=0.25)

res3 <- saltr2(alpha=0.05,data=darfur1,outcome=Y,treatment=D,bnch_reg=X,other_reg=X_oth,N=500)

}
