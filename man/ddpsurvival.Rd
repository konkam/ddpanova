\name{ddpsurvival}
\alias{ddpsurvival}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{DDP SURVIVAL}
\description{
  Initialize and run posterior MCMC for DDP ANOVA with right censored
  survival data.
}
\usage{
ddpsurvival(Y = NULL, D = NULL,
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
  All arguments are the same as in \code{\link{ddpanova}}, except the
  following 
  \item{Y}{(\code{n} by 3) data matrix (or file name).
    The first column is the event time. For rows with censored
    observations (2nd column=0), the first column is an initial value
    for an imputed event time.
    The second column is an indicator for an event, i.e., 0 for
    censored observations and 1 for event times.
    The third column is the last follow-up time
    (equal the 1st column for events).}
}
\details{
  See \code{\link{ddpanova-package}} for an explaination of
  all parameters.
  The function \code{ddpsurvival} initializes and carries
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
    ni <- 500
    y1   <-rnorm(ni, 3,.8)
    ## y2 = 0.6
    y21 <- rnorm(ni,2.0, 0.8)
    y22 <- rnorm(ni,4.0, 1.0)
    u <- runif(ni)
    y2 <- ifelse(u<0.6,y21,y22)
    y    <- c(y1,y2)
    nobs <- 2*ni

 ## now add censoring times
    yc <- rnorm(nobs,4,1.2)  
    event <- ifelse(y<yc,1,0)
    Y <- cbind(y,event,yc)

 ## design matrix with main effect mu and +/-offset
    d    <- cbind(rep(1,nobs),c(rep(-1,ni),rep(1,ni)))
    p    <- ncol(d)


  ## run MCMC
    d0 <- rbind(c(1,-1),c(1,1))  # design matrix for post predictive
    cat("\n Running 1000 iterations next -- be patient...\n")
    ddpsurvival(Y,d,n.iter=1000,d0=d0,ygrid=c(0,6))

  ## get posterior pred density for future observations
  ## with design matrix d=(1,-1) 
  ## rows of d0.
    pp <- post.pred()
    plot(pp$ygrid,pp$py[1,],type="l",xlab="Y",ylab="p",bty="l")
  ## overlay truth:
    p1 <- dnorm(pp$ygrid, 3.0, 0.8)
    lines(pp$ygrid, p1, lty=2)

  ## ... and d=(1,1) 
    plot(pp$ygrid,pp$py[2,],xlab="Y",ylab="p",bty="l",type="l")
    p2 <- 0.6*dnorm(pp$ygrid, 2.0, 0.8) +
          0.4*dnorm(pp$ygrid, 4.0, 1.0) 
    lines(pp$ygrid,p2, lty=2)

  ## and same as survival curves, assuming that the data were log(times)
    Sy  <- pp$Sy # survival curve
    S1  <- 1-pnorm(pp$ygrid, 3.0, 0.8)          # truth 1
    S2 <- 1 - (0.6*pnorm(pp$ygrid, 2.0, 0.8) +  # truth 2
               0.4*pnorm(pp$ygrid, 4.0, 1.0) )

    lygrid <- exp(pp$ygrid)                     # grid on time scale
    plot (lygrid,Sy[1,],type="l",bty="l",           # estimate
          xlab="Y",ylab="px",ylim=c(0,1),xlim=c(0,150))
    lines(lygrid,Sy[2,],type="l",col=2)             
    lines(lygrid,S1,    type="l",col=1,lty=2)       # simuulation truth
    lines(lygrid,S2,    type="l",col=2,lty=2)
}
}

