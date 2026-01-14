\name{ddpanova}
\alias{ddpanova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{DDP ANOVA}
\description{
  Initialize and run posterior MCMC for DDP ANOVA
}
\usage{
ddpanova(Y = NULL, D = NULL,
         n.iter = 1000, n.discard = 100, n.reinit = NULL,
         n.batch = 50, n.predupate = 100, n.printallpars = NULL,
         verbose = 3,
         m.prior = 1, B.prior = 1, S.prior=1, alpha.prior=2,
         s = NULL, S.init = NULL,
         q = NULL, R = NULL,
         m.init = NULL, B.init = NULL,
         cc = NULL, CC = NULL,
         a = NULL, A = NULL,
         alpha = 1, a0 = 1, b0 = 1,
         k0 = 1, member0 = 0,
         px = 1,
         d0 = NULL, resid = NULL,
         nx = 100, ygrid = NULL,
         xlist = 0, header = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{(\code{n} by 1) data vector (or file name)}
  \item{D}{(\code{n} by \code{p}) design matrix (or file name)}
  \item{n.iter}{number MCMC iterations}
  \item{n.discard}{initial burn-in}
  \item{n.reinit}{reinitialize every so many iterations.}
  \item{n.batch}{save posterior simulations every so many iterations}
  \item{n.predupate}{update posterior predictive every so many iterations}
  \item{n.printallpars}{print all parameters every so many iterations}
  \item{verbose}{0=silent}
  \item{m.prior}{indicator for resampling (vs. fixing) \eqn{m}.}
  \item{B.prior}{indicator for resampling (vs. fixing) \eqn{B}.}
  \item{S.prior}{indicator for resampling (vs. fixing) \eqn{S}.}
  \item{alpha.prior}{indicator for resampling (vs. fixing) \eqn{\alpha}.
    Use \code{alpha.prior=2} for fixing over the first
    \code{n.discard}/2 iterations, and imputing thereafter.}
  \item{s}{d.f. and scale in Gamma(\code{s}/2, 1/2*\code{s}*\code{S})
    base measure for precision (1/variance)}
  \item{S.init}{initial value for \code{S}}
  \item{q}{d.f. and scale in Gamma(\code{q},\code{q/R}) prior for \code{S}}
  \item{R}{(optional) prior scale parameter/other usage (see Details).}
  \item{m.init}{mean and variance in N(\code{m},\code{B}) base measure
    for location parameter}
  \item{B.init}{initial value for \code{B}}
  \item{cc}{degrees of freedom and matrix parameter for a \code{Wishart(cc, (1/cc) CC^{-1})} prior for \code{B^{-1}}.}
  \item{CC}{initial matrix for \code{CC}}
  \item{a}{mean and covariance matrix in N(\code{a},\code{A}) prior for
    \code{m}}
  \item{A}{initial matrix for \code{A}}
  \item{alpha}{initial value for the total mass parameter}
  \item{a0}{shape and scale of Gamma prior for \code{alpha}}
  \item{b0}{(optional) prior hyperparameter for \code{alpha}}
  \item{k0}{number of initial clusters}
  \item{member0}{(\code{n} by 1) vector of cluster membership
    indicators (only necessary if \code{k0} is not 1 or \code{n})}
  \item{px}{indicator for carrying out posterior predictive
    inference}
  \item{d0}{(\code{nd} times \code{p}) design matrix
    for \code{nd} future observations.}
  \item{resid}{binary vector of length \code{nd}, of indicators
    for including residual variance in the posterior predictive.}
  \item{nx}{grid size for plotting posterior predictive}
  \item{ygrid}{lower and upper bound of the grid for
    posterior predictive inference.}
  \item{xlist}{indicator for evaluating posterior fitted values
    for the observed data points.}
  \item{header}{indicator for the data files including
    a header line (only meaningful if \code{Y} and \code{D} are
    file names).}
}
\details{
  See \code{\link{ddpanova-package}} for an explaination of
  all parameters.
  The function \code{ddpanova} initializes and carries
  out MCMC posterior simulation. Simulation output is written to
  auxilary files (with .mdp extension) in the working directory.
  Use \code{\link{post.pred}} and \code{\link{post.sim}} to
  read in summaries.
}
\value{
  Posterior simulations are written into working files in the
  current working directory. Use \code{\link{setwd}} to change
  it if desired.

  The function has no return values.
}
\references{See \code{\link{ddpanova-package}}}
\author{Peter Mueller}
\examples{
\dontrun{
  ## simulate data from a mixture of two normal densities
    nobs <- 500
    y1   <-rnorm(nobs, 3,.8)
    ## y2 = 0.6
    y21 <- rnorm(nobs,2.0, 0.8)
    y22 <- rnorm(nobs,4.0, 1.0)
    u <- runif(nobs)
    y2 <- ifelse(u<0.6,y21,y22)
    y    <- c(y1,y2)

    ## design matrix with main effect mu and +/-offset
    d    <- cbind(rep(1,2*nobs),c(rep(-1,nobs),rep(1,nobs)))
    p    <- ncol(d)
    nobs <- 2*nobs


  ## run MCMC
    d0 <- rbind(c(1,-1),c(1,1))  # design matrix for post predictive
    ddpanova(y,d,n.iter=1000,d0=d0,ygrid=c(0,6))

  ## get posterior pred density for future observations
  ## with design matrix d=(1,-1)
  ## rows of d0.
    pp <- post.pred()
    plot(pp$ygrid[1,],pp$py[1,],type="l",xlab="Y",ylab="p",bty="l")
  ## overlay truth:
    p1 <- dnorm(pp$ygrid[1,], 3.0, 0.8)
    lines(pp$ygrid[1,],p1, lty=2)

  ## ... and d=(1,1)
    plot(pp$ygrid[2,],pp$py[2,],xlab="Y",ylab="p",bty="l",type="l")
    p2 <- 0.6*dnorm(pp$ygrid[2,], 2.0, 0.8) +
          0.4*dnorm(pp$ygrid[2,], 4.0, 1.0)
    lines(pp$ygrid[2,],p2, lty=2)
}
}

