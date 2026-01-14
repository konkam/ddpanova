/* ============================================================= 


/* ============================================================= 
/* ***********************************************************
   all MDP pars
/* *********************************************************** */

  struct MDP_PARS
{
  /* -----------------------------------------------------------
     Data & Parameters */
  /* ----------------------------------------------------------- */
  int censoring;                    /* indicator for censoring */
  double *Y,                        /* data */
    *CY,                            /* censoring time (if delta=0) */
    **d;                            /* (n x p) of design vec's  */
  int n_obs,                        /* n obs */
    *delta,                         /* indicator for censoring 
				       (0=cens;1=obs) */
    p;                              /* dimension of each obs vector */
  double *V, **mu,                /* V[i] and mu[i] */
    **mj,                           /* posterior mean for mu[j] */
    ***Vj;                          /* posterior cov matrix for mu[j] */
  int *updated;                     /* indicator if mj[j] is updated */
  double *V0;
  int  s;                           /* s, **S:  base distr of DP  */
  double S, R;
  int q;                            /* df of Wishhart hyperprior on S */
  double tau,                         
  w,                                /* prior shape for tau */
  W;                                /* prior scale for atu */
  double eta,                       /* latent par for sampling of alpha */
  alpha,                            /* conc par of DP */
  a0,                               /* prior shape for alpha */
  b0;                               /* prior scale for alpha */
  double  *mean,                      
  *aa, **A_inv;                     /* hyperpars on mean */
  double **B, **B_inv, **B_cd_inv,
  **C;
  int c;
  
  /* -----------------------------------------------------------
     indicators for expectations & predicitves to be computed *\
  /* ----------------------------------------------------------- */
       int px, pxy, mupxy, exy, exvecy, ex;       /* 0: no, 1: yes
					    the specific pars etc would then
					    be on files "pxy-init.mdp" etc. */
       int verbose; /* 1: full diagnostics, 0: not */
  /* -----------------------------------------------------------
     configuration */
  /* ----------------------------------------------------------- 
     The current configuration, i.e. number of distinct pi[i]'s 
     and classes of identical pi[i]'s is described by count, member
     and n_class.
     Notation: pi*[0],... pi*[n_class] is the list of distinct pi[i]'s.
     Of course n_class <= n_obs.
     n_class:                 number of distinct pi[i]'s.
     count[j]                 number of pi[i]'s equal to pi*[j]
     member[i]:               = j if pi[i] = pi*[j] 
     new[j]                   = r if class j was formed during
     iteration r (for debugging purposes) 
     first[j]                 = index of first obs in class j
     next[i]:                 = index of next obs in same class as obs i
     prev[i]:                 = index of previous obs in same class
     (prev[first] and next[last] = -1) 
     */ 
  int n_class, kmax,
  *count, *first,             
  *member, *next, *prev,
  *new;
};

/* -----------------------------------------------------------
   dummy function declarations
/* ----------------------------------------------------------- */
void  init();
void init_output();
void init_config();
void reinitialize_config();
void gibbs();
void finish();
void sample_config();
void check_class();
void take_out(int ind);
void add_class(int ind, int k);
void make_yhx(int k, double *y, double *h, double **x);
void make_Sk(int k, double *S);
void sample_mu(int k);
void sample_newV(); 
void sample_V(int k);
void new_class(int ind, double **V, double **V_cd_inv);
void sample_S();
void sample_m();
void sample_B();
void sample_alpha();
void print_pi();
void print_pars(int time);
int print_allpars();
void write_pars(int time);
int update_nclass();
int write_nclass();
void setup_mu(int k, int draw);

double loglik(int i, int l, double *m, double *s);
