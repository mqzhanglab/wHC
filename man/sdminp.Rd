% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdminp_calculation.R
\name{sdminp}
\alias{sdminp}
\title{This function calculates the Step Down minP method.}
\usage{
sdminp(res_ob, res_perm)
}
\arguments{
\item{res_ob}{is a numeric vectors of statistics}

\item{res_perm}{is a matrix with each colum is a permuated statistics with the same length as res_ob}
}
\value{
a numeric vector adjusted pvalues
}
\description{
It implement the improved step-down minP algorithm by Ge et al, 2003 (box 4)
}
\examples{
ob=rnorm(10,2,2)
perm=matrix(rnorm(100,2,2),ncol=10,nrow=10)
sdminp(ob,perm)
}
\references{
Ge, Y., Dudoit, S., & Speed, T. P. (2003). Resampling-based multiple testing for microarray data analysis. Test, 12(1), 1–77. https://doi.org/10.1007/BF02595811
}
