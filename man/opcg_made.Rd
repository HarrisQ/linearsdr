% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/opcg_wrap_cpp.R
\name{opcg_made}
\alias{opcg_made}
\title{OPCG-MADE - Local gradient estimation}
\usage{
opcg_made(
  x_matrix,
  y_matrix,
  bw,
  lambda,
  B_mat = NULL,
  ytype = "continuous",
  method = "newton",
  parallelize = F,
  r_mat = NULL,
  control_list = list()
)
}
\arguments{
\item{x_matrix}{a 'nxp' matrix of predictors;}

\item{y_matrix}{a 'nxm' response;}

\item{bw}{the bandwidth parameter for the kernel; the default kernel is gaussian}

\item{lambda}{an L2 penalty term for the negative log-likelihood}

\item{B_mat}{the fixed coefficient matrix in MADE-step of MADE;
not needed for OPCG, i.e. is set to the identity}

\item{ytype}{the response type; continuous, categorical or ordinal}

\item{method}{"newton" or "cg" methods; for carrying out the optimization using
the standard newton-raphson (i.e. Fisher Scoring) or using Congugate Gradients}

\item{parallelize}{Default is False; to run in parallel, you will need to have
foreach and some parallel backend loaded; parallelization is strongly recommended
and encouraged.}

\item{r_mat}{a 'pxd' matrix for refining the weights in rOPCG and rMADE}

\item{control_list}{a list of control parameters for the Newton-Raphson
or Conjugate Gradient methods}
}
\value{
\itemize{
\item ahat - List of estimated local intercepts
\item Dhat - List of estimated local slopes/gradients
\item Dhat_ls - List of initial values for local slopes/gradients;
for least squares, these are the same as the Dhat
\item weights - The kernel weights used in the local-linear estimation;
}
}
\description{
This is an internal function called by OPCG. MADE also uses this
function in its OPCG-step. This estimates the local intercept and
slope coefficients.
}
