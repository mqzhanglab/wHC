% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/higher_criticism.R
\name{cal_cdf}
\alias{cal_cdf}
\title{This function calculates the weighted-pvalues.}
\usage{
cal_cdf(pm, w = 1)
}
\arguments{
\item{pm}{is a statistics matrix of P-values from nrow=n genes(independent tests),ncol=d. Pm are not encouraged to have only 1 rows, if that happend, warning massage will produced and this function returns sqrt(1/pm-1), rather than the double side statistics from linear regression.Thus, we have S(t)~Binomial(n,p).}

\item{w}{is the weight, if not specify, w=1, if specify w must have the same length as nrow Pm}
}
\value{
numeric vector with each elements is a Higher criticism values calculated from each colum of the Pm
}
\description{
It calculates the CDF(Cumulative Distribution Function) of weighted pvalues from the original pvalues and their weights and assign new pvalues based on this CDF.
}
\examples{
pval=matrix(runif(200,0,1),ncol=4,nrow=50)
w0=seq(0.5,1.5,length=50)
pwval=cal_cdf(pval,w=w0)
}
\references{
Genovese, C. R., Roeder, K., & Wasserman, L. (2006). False Discovery Control with p-Value Weighting. Biometrika, 93(3), 509–524.
}
