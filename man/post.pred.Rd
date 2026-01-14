\name{post.pred}
\alias{post.pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Posterior predictive inference }
\description{
  Return posterior predictive densities for future observations.
}
\usage{
post.pred()
}
%- maybe also 'usage' for other objects documented here.
\details{
  Must have called \code{\link{ddpanova}} or
  \code{\link{ddpsurvival}} first.
  
  Returns the posterior predictive distribution
  \eqn{p(y^* \mid d^*, data)} for a future observation
  \eqn{y^*=y_{n+1}} with design vector \eqn{d^*}.
  The argument \code{d0} of \code{\link{ddpanova}} specifies the
  list of \code{nd} future designs that are considered.
}
\value{
  The function returns a list:
  \item{ygrid}{(\code{nx} x 1) vector 
    of grid values for the predictive
    for future units with design vectors given in \code{d0}}.
  \item{py}{(\code{nd} by \code{nx}) matrix of
    posterior predictive probabilities \eqn{p(y^* \mid d^*,data)},
    evaluated on the grid given in \code{ygrid}.}
  \item{Sy}{(\code{nd} by \code{nx}) matrix of
    estimated survival curves for the \code{nd} future units with
    design vectors specified in the rows of \code{d0}.
    Only when called after \code{\link{ddpsurvival}}.}
  \item{sSy}{(\code{nd} by \code{nx}) matrix of
    pointwise posterior predictive s.d.'s for the survival curve.
    Only when called after \code{\link{ddpsurvival}}.}
  \item{hy}{(\code{nd} by \code{nx}) matrix of
    posterior predictive hazards for future units with design
    vectors specified in \code{d0}.
    Only when called after \code{\link{ddpsurvival}}.}
  \item{shy}{(\code{nd} by \code{nx}) matrix of pointwise s.d.'s of
    the hazard curves.}
  \item{censoring}{indicator whether call was following
    \code{\link{ddpsurvival}} (\code{censoring=1}) or
    \code{\link{ddpanova}} (\code{censoring=0}).}
}
\references{see \code{\link{ddpanova-package}}}
\author{Peter Mueller}
\note{Use \code{\link{post.sim}} to obtain a Monte Carlo sample
  of parameter values.

  See \code{\link{ddpanova}} for an example.
}
