double lpost(double *theta, struct Model mdl, struct Data dta);


/* *************************************************
   Model structure 
 ************************************************* */
struct Model{
  int 
    lag_mean,      /* lags in mean equation, =1 currently */
    lag_var;
  int 
    p,             /* total n pars */
    pb,pg,         /* length of beta and gamma vector */
    pg_lags, pg_crosslags, pg_wkdays;
                   /* length of compontents of gamma */
  int
    p1;            /* number of pars to be printed each iterations */
};

/* *************************************************
   Data structure
 ************************************************* */
struct Data {
  int d;             /* dimension of data vector at time t */
  int n;             /* number observations */
  double **y;        /* data matrix */
};

/* *************************************************
   Markov chain Monte Carlo pars
 ************************************************* */
struct MCMC {
  int iter,         /* number of the current iteration */
  niter,            /* number iterations */
  nbatch,           /* batch siyes for each step of indep chain */
  npi,              /* size of batches for ergodic avges */
  ndiscard,         /* initial burn-in */
  ncollect,         /* start accumulating MC sample after ncollect its */
  nrepar,           /* reparametrize every nrepar its */
  nstoprepar,       /* stop reparametrization after nstoprepar its */
  nlastrepar,       /* time of last reparametrization */
  nmetr,            /* run metr every nmetr iterations */
  nbatchmetr,       /* size of "mini metr" */
  verbose,          /* be verbose about comments */
  seed,             /* seed for RV generator */
  t0;               /* par for each reparametrization */
  double
    c_indep,
    c;              /* scale factor for Metropolis */
};

/* *************************************************
   Monte Carlo sample
 ************************************************* */
struct MCsample {
  double 
    **theta;       /* Monte Carlo sample, each row is one par vector */
  int
    m;
};

/* *************************************************
   Reparametrization
 ************************************************* */
struct Repar {
  double *mu, *sd;
};

/* *************************************************
   Erogidc averages
 ************************************************* */
struct Ergodic {
  double	 
    *th,       /* ergodic mean of th's */
    *th2,      /* ergodic mean of th^2 */
    *sd,       /* estimated posterior sd's */
    *thbar,    /* mean over parallel m chains at current iteration */
    *sdbar;    /* sd estimated over m parallel chains */
  double
    nupdate;  /* number of updates for ergodic avgs */
};

  
