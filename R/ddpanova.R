
'ddpsurvival' <- function(Y=NULL,  # responses, (n x 1) vector
                                ## or (n x 3) for censored (or file name)
                     D=NULL,    # design matrix, (n x p) matrix (or file nm)
                     censoring=1,    # indicator for censoring
                     n.iter=1000,    # n MCMC iterations
                     n.discard=100,  # initial burn in
                     n.reinit=NULL,  # re-inialize every so many iterations
                     n.batch=25,     # save posterior sim's every so many
                     n.predupdate=25,# update posterior pred sim ..
                     n.printallpars=NULL, # save all imputed pars ..
                     ny0.print=25,   # save trajectories for first 10 imputed censored times
                     verbose=3,      # 0=silent
                     m.prior=1,   # indicator for resampling m
                     B.prior=1,   # .. and B
                     S.prior=0,   # .. and S
                     alpha.prior=2,   # .. and S
                     s=NULL,      # d.f. for basemeasure G0(1/V)=Ga(s/2,sS/2)
                     S.init=NULL, # initial value for S
                     q=NULL,      # d.f. for prior S ~ Ga(q,q/R)
                     R=NULL,      # ..
                     m.init=NULL, # initial values for G0(mu)=N(m,B)
                     B.init=NULL, # ..
                     cc=NULL,     # d.f. for prior 1/B ~ Wish(c,1/(cC))
                     CC=NULL,     # matrix par..
                     a=NULL,      # prior moments for m~(a,A)
                     A=NULL,      # ..
                     alpha=1,     # initial value for total mass par
                     a0=1,        # prior pars for alpha~Ga(a0,b0)
                     b0=1,
                     k0=NULL,      # initial number of clusters
                     member0=0,    # initial cluster indicators
                                   ## if k0!=1 and k0!=n
                     ## parameters for posterior predictive
                     px=1,           # indicator for computing
                                     ## posterior predictive inference
                     d0=NULL,     # (nd x p) matrix of designs for
                                  ## posterior pred inference
                     nd0=0,       ## save simulations for the first nd0
                                  ## future designs in d0
                     resid=NULL,  # (nd x 1) vector of indicators
                                  ## to include residual variance in pred
                     nx=100,      # grid size for posterior pred
                     ygrid=NULL,  # (nd x 2) matrix of lo & hi for
                                  ## posterior predictiv grids
                     xlist=0,     # indicator for posterior pred
                                  ## inference for data designs
                                  ## of observed data points (w/o resid)
                     sim=0,       ## posterior pred simulations
                     nsim=10,
                     header=0)    # indicator for whether data
                                  ## files Y and D include headers

{
  if ( (censoring==1) & (!is.matrix(D)) ){
    cat("\n *** Error: when using censoring (event time data),\n")
    cat("\n     use data matrix with 3 columns. See help.\n")
    return(-1)
  }
  return(
         ddpanova(Y,D,censoring,
                  n.iter,n.discard,n.reinit,n.batch,n.predupdate,
                  n.printallpars,ny0.print,verbose,
                  m.prior,B.prior,S.prior,alpha.prior,
                  s,S.init, q,R, m.init, B.init, cc, CC,
                  a,A,alpha,a0,b0,k0,member0,
                  px,d0,nd0,resid,nx,ygrid,xlist,sim,nsim,header)
         )
}

  
'ddpanova' <- function(Y=NULL,  # responses, (n x 1) vector
                                ## or (n x 3) for censored (or file name)
                     D=NULL,    # design matrix, (n x p) matrix (or file nm)
                     censoring=0,    # indicator for censoring
                     n.iter=1000,    # n MCMC iterations
                     n.discard=100,  # initial burn in
                     n.reinit=NULL,  # re-inialize every so many iterations
                     n.batch=25,     # save posterior sim's every so many
                     n.predupdate=25,# update posterior pred sim ..
                     n.printallpars=NULL, # save all imputed pars ..
                     ny0.print=0,   # save trajectories for first 10 imputed censored times
                     verbose=3,      # 0=silent
                     m.prior=1,   # indicator for resampling m
                     B.prior=1,   # .. and B
                     S.prior=1,   # .. and S
                     alpha.prior=2,   # .. and S
                     s=NULL,      # d.f. for basemeasure G0(1/V)=Ga(s/2,sS/2)
                     S.init=NULL, # initial value for S
                     q=NULL,      # d.f. for prior S ~ Ga(q,q/R)
                     R=NULL,      # ..
                     m.init=NULL, # initial values for G0(mu)=N(m,B)
                     B.init=NULL, # ..
                     cc=NULL,     # d.f. for prior 1/B ~ Wish(c,1/(cC))
                     CC=NULL,     # matrix par..
                     a=NULL,      # prior moments for m~(a,A)
                     A=NULL,      # ..
                     alpha=1,     # initial value for total mass par
                     a0=1,        # prior pars for alpha~Ga(a0,b0)
                     b0=1,
                     k0=NULL,      # initial number of clusters
                     member0=0,    # initial cluster indicators
                                   ## if k0!=1 and k0!=n
                     ## parameters for posterior predictive
                     px=1,           # indicator for computing
                                     ## posterior predictive inference
                     d0=NULL,     # (nd x p) matrix of designs for
                                  ## posterior pred inference
                     nd0=0,       ## save simulations for the first nd0
                                  ## future designs in d0
                     resid=NULL,  # (nd x 1) vector of indicators
                                  ## to include residual variance in pred
                     nx=100,      # grid size for posterior pred
                     ygrid=NULL,  # (nd x 2) matrix of lo & hi for
                                  ## posterior predictiv grids
                     xlist=0,     # indicator for posterior pred
                                  ## inference for data designs
                                  ## of observed data points (w/o resid)
                     sim=0,       ## posterior pred simulations
                     nsim=10,
                     header=0)    # indicator for whether data
                                  ## files Y and D include headers

{
  ## read in data
  options(digits=2)
  
  if(is.null(Y) | is.null(D)){
    cat("\n *** Error: need data vector (or file name).\n")
    return(-1)
  }
  if(!is.element(header,c(0,1))){
    cat("\n *** Error: header indicator must be 0 or 1.\n")
    return(-1)
  }
  if(!is.element(censoring,c(0,1))){
    cat("\n *** Error: censoring indicator must be 0 (no) or 1 (censoring).\n")
    return(-1)
  }
  if (is.numeric(Y)){ # input is response values already
    y <- Y
  } else { # read from data file
    if (censoring==0){
      y <- scan(Y,skip=header)
    } else {
      y <- as.matrix(read.table(Y,skip=header))
    }
  }
  if (censoring==0){
    ny <- length(y)
    yvec <- y
  }  else {
    ny <- nrow(y)
    yvec <- y[,1]
  }
  if (is.numeric(D)){ # input is response values already
    d <- D
  } else { # read from data file
    d <- as.matrix(read.table(D,header=header))
  }
  if(ny != nrow(d)){
    cat("\n *** Error: nrow(design) != length(y).\n")
    cat("\n ny, nrow(d)=",ny,nrow(d),"\n")
    return(-1)
  }
  n <- nrow(d)
  p <- ncol(d)
  ## MCMC 
  if (n.discard >= n.iter){
    cat("\n *** Error: n.iter <= n.discard. \n")
    return(-1)
  }
  if (n.batch >= n.iter/2){
    cat("\n *** Error: n.iter/2 <= n.batch. \n")
    return(-1)
  }
  if (!is.element(verbose,0:3)){
    cat("\n *** Error: verbose not in {0,1,2,3}. \n")
    return(-1)
  }
  if(is.null(n.reinit))
    n.reinit <- n.iter+1
  if(is.null(n.printallpars))
    n.printallpars <- n.iter+1
  
  ## initial values for pars and hyperparameters
  XtX <- t(d)%*%d
  if(min(eigen(XtX)$values) < 0.001){ # colinearity
    cat("\n *** Error: design matrix almost colinear.\n")
    return(-1)
  }
  H <- solve(XtX)       # hat matrix for regression..
  bhat <- H%*%t(d)%*%yvec  # mle of single regression - could be aweful!
  yhat <- d%*%bhat      # fitted values
  e <-  yvec-yhat          # residuals 
  sig2 <- c(var(e))     # residual varance
  Sigma <- sig2*H       # mle covariance matrix
  S0 <- 4*diag( round10(diag(Sigma)) )
  ## base measure G0(1/V[i]) = Ga(s/2, 1/2*s*S)
  if(!is.null(S.init)){
    if(length(S.init)!=1){
      cat("\n *** Error: S.init must be a scalar.\n")
      return(-1)
    }
  } else 
    S.init <- sig2
  if(is.null(s)){
    s <- 4
  }
  ## base measure G0(mu[i]) = N(m,B) 
  if(!is.null(B.init)){
    if(any(dim(B.init) != c(p,p))){
      cat("\n *** Error: B.init must be a (p x p) matrix.")
      cat("\n     p=ncol(design matrix).\n")
      return(-1)
    }
  } else 
    B.init <- 4*Sigma
  if(!is.null(m.init)){
    if(length(m.init) != p){
      cat("\n *** Error: dim(m.init) != p\n.")
      return(-1)
    }
  } else 
    m.init <- bhat
  ## hyperprior 1/B ~ Wish(c, 1/c*C^-1)
  if(!is.null(CC)){
    if(any(dim(CC) != c(p,p))){
      cat("\n *** Error: C must be a (p x p) matrix.")
      cat("\n     p=ncol(design matrix).\n")
      return(-1)
    }
  } else 
    CC <- S0
  if(is.null(cc))
    cc <- p+5
  ## hyperprior on m; m~N(a,A)
  if(!is.null(A)){
    if(any(dim(A) != c(p,p))){
      cat("\n *** Error: A must be a (p x p) matrix.")
      cat("\n     p=ncol(design matrix).\n")
      return(-1)
    }
  } else 
    A <- S0
  A.inv <- solve(A)
  if(!is.null(a)){
    if(length(a) != p){
      cat("\n *** Error: dim(a) != p\n.")
      return(-1)
    }
  } else 
    a <- bhat
  ## hyperprior on S, S ~ Ga(q,q/R)
  if(!is.null(R)){
    if(length(R)!=1){
      cat("\n *** Error: R must be a scalar.\n")
      return(-1)
    }
  } else 
    R <- round10(sig2)
  if(is.null(q))
    q <- 4
  ## initial cluster membership
  if (is.null(k0)){ # 3 default clusters for low, medium & high resids
    smy <- quantile(yvec,c(0.25,0.75))
    member0 <- ifelse(yvec<smy[1],0,
                      ifelse(yvec<smy[2],1,2))
    k0 <- 3
    ## member0 <- 0
    ## k0 <- 1
  }
  if ((k0!=1) & (k0!=n) & (is.null(member0))){
    cat("\n *** Error: when specifying 1 < k0 < n,")
    cat("\n            then need initial cluster membership member0.\n")
    return(-1)
  }
  if ((k0!=1) & (k0!=n) & (any(!is.element(member0,0:(k0-1))))){
    cat("\n *** error: member0 needs to be in {0...k0-1}.\n")
    return(-1)
  }
  
  initfile <- "init.mdp"
  cat(" n ",ny,"\n",
      " p ",p,   "\n",
      " censoring ", censoring, "\n",
      " niter ", n.iter, "\n",
      " m_prior ", m.prior, "\n",
      " B_prior ", B.prior, "\n",
      " S_prior ", S.prior, "\n",
      " alpha_prior ", alpha.prior, "\n",
      " ndiscard ", n.discard, "\n",
      " nreinit ", n.reinit, "\n",
      " npi ", n.batch, "\n",
      " npredupdate ", n.predupdate, "\n",
      " nprintallpars ", n.printallpars, "\n",
      " ny0 ", ny0.print, "\n",
      " verbose ", verbose, "\n",
      " px ", px, "\n",
      " s ", s, "\n",
      " S-init ", S.init, "\n",
      " q ", q, "\n",
      " R ", R, "\n",
      " B-init \n", rbind(format(t(B.init)),"\n"),
      " c ", cc, "\n",
      " C ", rbind(format(t(CC)),"\n"), "\n",
      " m-init \n ", format(m.init),"\n",
      " aa \n ", format(a), "\n",
      " A-inv \n ", rbind(format(t(A.inv)),"\n"),
      " alpha-init ", alpha, "\n",
      " a0 ", a0, "\n",
      " b0 ", b0, "\n",
      " k0 ", k0, "\n",
      " member0 \n", member0, "\n",
      file=initfile)
  ## write data matrices
  dtafile <- "data.mdp"
  write(format(t(y)),ncolumns=3,file=dtafile)
  dtafile <- "data-d.mdp"
  write(format(t(d)),ncolumns=p,file=dtafile)

  
  ## write px-init.mdp
  if(px==1){
    if (is.null(d0)){
      cat("\n *** Error: need to give (nd x p) design matrix d0")
      cat("\n            for posterior predictive or set px=0.\n")
      return(-1)
    }
    if(ncol(d0) != p){
      cat(dim(d0),p,"\n")
      cat("\n *** Error: ncol(d0) != p.\n")
      return(-1)
    }
    nd <- nrow(d0)
    if (nd0 > nd){
      cat("\n *** Error: need nd0 <= nd.\n")
      return(-1)
    }
    if (is.null(resid))
      resid <- rep(1,nd)
    if (length(resid)==1)
      resid <- rep(resid,nd)
    if (length(resid)!=nd){
      cat("\n *** Error: length(resid) != nrow(d0).\n")
      return(-1)
    }
    if (nx<25){
      cat("\n *** Error: grid size nx should be at least 25.\n")
      return(-1)
    }
    if (is.null(ygrid))
      ygrid <- range(yvec)
    if (length(ygrid) != 2){
      cat("\n *** Error: ylist needs to be a pair of min and max for the grid.\n")
      return(-1)
    }
    if(!is.element(xlist,0:1)){
      cat("\n *** Error: xlist must be 0 or 1.\n")
      return(-1)
    }
    if(!is.element(header,0:1)){
      cat("\n *** Error: header must be 0 or 1.\n")
      return(-1)
    }
    initfile <- "px-init.mdp"
    cat("  nx ",nx,"\n",
        " nd ",nd,"\n",
        " nd0 ",nd0,"\n",
        " sim ", sim, "\n",
        " nsim ", nsim, "\n",
        " d \n",rbind("  ",format(t(d0)),"\n"),
        " resid ", format(resid), "\n",
        " ygrid ", format(ygrid),"\n",
        " xlist ", 0,
        file=initfile)
  } ## px==1
  .C("ddpanova",package="ddpanova")
}



'post.pred' <- function()
{ ## returns grid and estimated post pred density for post predictive
  
  xgridfile <- "px-x.mdp"
  pfile <- "px-p.mdp"
  Sfile <- "px-Sy.mdp"
  ## geet nx and nd from 'px-init.mdp'
  inp <- scan(file="px-init.mdp",nmax=2,what=list("nx",1),
              comment.char="#")
  nx <- inp[[2]][1]
  nd <- inp[[2]][2]

  ygrid <- scan(xgridfile)
  px <- matrix(scan(pfile),ncol=nx,byrow=T)
  if(file.exists(Sfile)){
    censoring <- 1
    Sx <- matrix(scan(Sfile),ncol=nx,byrow=T)
    sSx <- matrix(scan("px-sSy.mdp"),ncol=nx,byrow=T)
    hx <- matrix(scan("px-hy.mdp"),ncol=nx,byrow=T)
    shx <- matrix(scan("px-shy.mdp"),ncol=nx,byrow=T)
  } else {
    censoring <- 0
    Sx <- NULL
    sSx <- NULL
    hx <- NULL
    shx <- NULL
  }
  
  if(length(ygrid) != nx){
    cat("\n *** Error: length(px-x.mdp) != nx.")
    cat("\n            Make sure you run MCMC simulation first in")
    cat("\n            the same directory, using 'ddpanova()'.\n")
    return(-1)
  }
  if(any(dim(px) != c(nd,nx))){
    cat("\n *** Error: dim(px-p.mdp) != (nd x nx).")
    cat("\n            Make sure you run MCMC simulation first in")
    cat("\n            the same directory, using 'ddpanova()'.\n")
    return(-1)
  }
  return(list(ygrid=ygrid,py=px,Sy=Sx,sSy=sSx,hy=hx,shy=shx,
              censoring=censoring))
}

'post.sim' <- function()
{ # return posterior simuulations for anova parameters
  zfile <- "z.mdp"
  ## geet n and p 'init.mdp'
  inp <- scan(file="init.mdp",nmax=2,what=list("n",1),
              comment.char="#")
  n <- inp[[2]][1]
  p <- inp[[2]][2]
  Z <- matrix(scan(zfile),byrow=T,ncol=p)
  return(Z)
}

round10 <- function(x)
  {
    return(exp(log(10)*ceiling(log10(x))))
  }
    

## used for debugging only...
##
initial.values <- function()
{
  d0 <- rbind(c(1,-1),c(1,1))  # for prediction
  censoring <- 0
  header <- 0; # Y <- "data.mdp"; D <- "data-d.mdp"
  m.prior <- 1; B.prior <- 1; n.iter <- 500; n.discard <- 100
  n.reinit <- n.iter
  n.batch <- 50; n.predupdate <- 100; n.printallpars <- n.iter
  verbose <- 3; px <- 1;
  m.init <- NULL; B.init <- NULL
  s <- 5; S.init <- NULL
  q <- 5; R <- NULL
  cc <- 5; CC <- NULL
  A <- NULL; a <- NULL
  alpha <- 1; a0 <- 1; b0 <- 1
  nx <- 50; nd <- 2; sim <- 1; nsim <- 10;
  d <- rbind(c(1,-1),c(1,1))
  resid <- c(1,1)
  ygrid <- c(0,6); xlist <- 0
  S.init <- NULL
  k0 <- 1; member0 <- 0
}  
