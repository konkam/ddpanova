\name{post.sim}
\alias{post.sim}
\title{Posterior Monte Carlo sample}
\description{
  Retrieves a simulated samples \eqn{\mu} from the posterior estiamted
  mixing measure \eqn{\bar G = E(G \mid data)}.
}
\usage{
  post.sim()
}
%- maybe also 'usage' for other objects documented here.
\details{
  The function returns a matrix of p-dimensional vectors
  \eqn{\mu \sim \bar G}. See \code{\link{ddpanova-package} for a
    statement of the model and the parameters \eqn{\mu}.
    Alternatively, one can think of \eqn{\mu} as posterior predictive
    sample \eqn{\mu \sim \int G(\mu) dp(G \mid data)}.
    The samples can be used to evaluate arbitrary posterior and/or
    posterior predictive summaries.
}
\value{
  The function returns a matrix
  \item{Z }{(\code{M} by \code{p}) matrix with one simuulation in each line.}
}
\references{see \code{\link{ddpanova-package}}}
\author{Peter Mueller}
\note{see  \code{\link{ddpanova}} for an example.}
