% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr.genius.R
\name{data_gen_x}
\alias{data_gen_x}
\title{Simulate a dataset with covariates}
\usage{
data_gen_x(m, n, beta0, gamma)
}
\arguments{
\item{m}{Number of SNPs}

\item{n}{Sample size}

\item{beta0}{True causal effect}

\item{gamma}{A parameter that controls the magnitude of heteroscedasticity, i.e., the identification strength}
}
\value{
A list
\describe{
\item{z}{A n*m matrix of SNPs.}
\item{a}{A n-dimentional vector for the exposure.}
\item{y}{A n-dimentional vector for the outcome.}
\item{x}{A n*2 matrix for the covariate.}
}
}
\description{
Simulate a dataset with covariates
}
\examples{
df<-data_gen_x(m=100,n=1e5,beta=0.4,gamma=0.1)
}
\references{
Ting Ye, Zhonghua Liu, Baoluo Sun, and Eric Tchetgen Tchetgen (2021). GENIUS-MAWII: For Robust Mendelian Randomization with Many Weak Invalid Instruments.\url{https://arxiv.org/abs/2107.06238}.
}
