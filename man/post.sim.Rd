\name{post.sim}
\alias{post.sim}
\title{Posterior Monte Carlo sample}
\description{
  Retrieves simulated samples \eqn{\mu} from the posterior-estimated
  mixing measure \eqn{\bar G = E(G \mid data)}.
}
\usage{
  post.sim()
}
%- maybe also 'usage' for other objects documented here.
\details{
  The function returns a matrix of p-dimensional vectors
  \eqn{\mu \sim \bar G}. See \link{ddpanova-package} for a
  statement of the model and the parameters \eqn{\mu}.
  Alternatively, one can think of \eqn{\mu} as a posterior predictive
  sample: \eqn{\mu \sim \int G(\mu)\, dP(G \mid data)}.
  The samples can be used to evaluate arbitrary posterior and/or
  posterior predictive summaries.
}
\value{
  A matrix \item{Z}{An \code{M} by \code{p} matrix with one simulation in each row.}
}
\references{see \code{\link{ddpanova-package}}}
\author{Peter Mueller}
\note{see  \code{\link{ddpanova}} for an example.}
