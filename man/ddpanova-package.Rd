\name{ddpanova-package}
\alias{ddpanova-package}
\docType{package}
\title{
DDP ANOVA with univariate response
}
\description{
  The package implements posterior simulation and
  posterior predictive inference for the DDP ANOVA model, assuming
  a univariate response.
}
\details{
  \tabular{ll}{
    Package: \tab ddpanova\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2007-04-23\cr
    License: \tab GNU?\cr
  }
  Use \code{\link{ddpanova}} to initialize and run the posterior
  MCMC simulation. Use \code{\link{post.pred}} to get posterior
  predictive densities for the observable outcomes under assumed
  design vectors for future patients.
  Use \code{\link{post.sim}} to get posterior simulations for the ANOVA
  parameters.

  The model is as described in De Iorio et
  al. (2004).
  The same model is used in Mueller et al. (2005) as the random effects
  distribution in a repeated measurements model.
  Let \eqn{y_i} denote the i-th observed response for an individual
  (experimental unit) with design vector \eqn{d_i}.
  For example, if the model is used to jointly analyze data from two
  related clinical trials, with 3 treatment levels in each trial,
  the design vector for a patient in study 2, at treatment level 3 might
  be \eqn{d_i=(1,-1,0,1)}, including an intercept, an offset for the 2nd
  study, no offset for treatment level 2 and an offset for treatment
  level 3.

  Let \eqn{N(x; m,s), Ga(x; a,b)} and  \eqn{Wish(f,A)}
  denote a normal p.d.f. with moments \eqn{(m,s)}, a Gamma p.d.f. with
  shape \eqn{a} and expectation \eqn{a/b}, and a Wishart distribution with
  degrees of freedom \eqn{f} and expectation \eqn{f A}.

  We assume a mixture of ANOVA sampling model
  \deqn{
    y_i \sim \int N(y_i; d_i^\prime\mu, V) dG(\mu,V),
  }
  with a Dirichlet process (DP) prior for the random mixing measure,
  \eqn{G \sim DP(G_0,\alpha)}. Here, \eqn{G_0(\mu,V)} is the base measure and
  \eqn{\alpha} is the total mass parameter.

  Let \eqn{H_d(y) = \int N(y; d^\prime\mu,V) dG(\mu,V)} denote the sampling model
  for an observation with design vector \eqn{d}. The ANOVA DDP defines a
  joint probability model on the set of random probability measures
  \deqn{
    \{H_d(\cdot),\; d \in D\}
  }
  for designs from a design space \eqn{D}.
  The model is completed with conditioally conjugate hyperpriors.
  We assume \eqn{G_0(\mu,V) = Ga(1/V; s/2, sS/2) N(\mu; m,B)} and
  hyperpriors \eqn{S \sim Ga(0.5 q, 0.5 q/R), m \sim N(a,A)},
  \eqn{B^{-1} \sim Wish(c,(cC)^{-1})} and
  \eqn{\alpha \sim Ga(a_0,b_0)}.
  In summary
  \deqn{
    \begin{eqnarray*}
    y_i & \sim & \int N(y_i;\; d_i^\prime\mu, V)\; dG(\mu,V),\\
    G   & \sim & DP(\alpha,G_0),
    G_0(\mu,V)  = Ga(1/V; s/2, sS/2) N(\mu; m,B)\\
    S          & \sim & Ga(q,q/R),  E(S)=R\\
    m          & \sim & N(a,A)\\
    B^{-1}     & \sim & Wish(c,(cC)^{-1}),  E(B^{-1})=C^{-1}\\
    \alpha      & \sim & Ga(a_0,b_0).
    \end{eqnarray*}
  }
}
\author{
Peter Mueller

Maintainer: Peter Mueller <pm@wotan.mdacc.tmc.edu>
}
\references{
  De Iorio, M., Mueller, P., Rosner, G., and Maceachern, S. (2004).
  \emph{An ANOVA Model for Dependent Random Measures},
  Journal of the American Statistical Association, 99(465), 205--215.

  Mueller, P., Rosner, G., De Iorio, M., and MacEachern, S. (2005).
  A Nonparametric Bayesian Model for Inference
    in Related Studies.
  \emph{Applied Statistics}, 54 (3), 611-626.
}
\keyword{ package }
