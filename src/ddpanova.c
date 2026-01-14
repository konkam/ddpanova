/* ============================================================= 
   DDP ANOVA with univariate outcome
   
/* ============================================================= 
   Model:
   y[i] | pi[i] is N(y[i]; mu[i]'*d[i], V[i])
                   pi[i] = (mu[i], V[i])
   pi[i]        is DP(G0, alpha)
   G0(pi)       is Ga(1/V; s/2, 1/2*s*S) N(mu; m,B), 
                   i.e. E(1/V)=1/S
		        use s > p+3 (for 2nd moments)
   S            is Ga(S; q, q/R), i.e. E(S)=R
                p(S | V) = Ga(q+k*s; ....)
   B            W(1/B; cc; 1/cc*1/C); s.t. E(1/B)=1/C		
   alpah        is Ga(a0, b0), E(alpha) = a0/b0

Estimation as described in MacEachern & Peter (1994)
Notation:
   Some of the pi[i]'s can be identical. The list of distinct
   pi[i]'s is denoted by pi*[0]..pi*[n_class].

sanNoninformative hyperprior pars:
   S:     R_inv=0, q=p+4
   mean:  A_inv=0, aa=any
   alpha: a0=0, b0=0 

   Details as in Escobar & West */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "ddpanova.h"
#include "bayes.h"
#include "rand.h"
#include "nrutil.h"
#include "matrix.h"
#include "vector.h"
#include "rand-nr.h"
#include "mess.h"
#include "interface.h"
#include "predict.h"
#include "regr.h"

#define Ndiscard 50          /* number of discarded iterations */

/* ***********************************************************
   all MDP pars
/* *********************************************************** */
struct MDP_PARS mdp;
int 
  q=1, /* dimension of response; i.e., DIMENSION of ANOVA eff = r*q */
       /* NOTE: currently only set up for q=1 */
  p, /* dimension of ANOVA effects = r*q */
  r; /* NUMBER of ANOVA effects */  

double *m_aux, *s_aux; /* auxiliary arrays of dim nobs
			   for temporary results */

/* -----------------------------------------------------------
   Gibbs iterations
/* ----------------------------------------------------------- */
int n_iter,             /* max iterations */
  iter,                 /* current iter number */
  n_discard,            /* discard the first n_discard draws */
  n_reinitialize,         /* reinitialize configurations */
  n_predupdate,
  n_pi,                 /* every n_pi update expectations etc. */
  ny0_print;             /* number of imputed (censored) obs to save
			   for convergence diagnostics 
			   (only when censoring==1) */

int n_printallpars;     /* intervals to print all pars if 
			   verbose > 2 */
int m_prior, B_prior, S_prior, alpha_prior;
int n_class0, *member0; /* initial config */

/* -----------------------------------------------------------
   various output files
/* ----------------------------------------------------------- */
FILE *parFile, *muFile, *VFile, *SFile, *meanFile, *BFile, *y0file;

long is1=123456789, is2=981963; 

/* -----------------------------------------------------------
   debugging
/* ----------------------------------------------------------- */
int dbg=0;

/* ============================================================= 
   main
/* ============================================================= */ 
void ddpanova()  /* for calling as a function.. */
/* void main()  /* for debugging.. */
{ 
  init();
  init_output();
  pred_init(mdp);
  init_config();
  gibbs(); 
  finish();
}				  


/* ============================================================= 
   functions called from main
/* ============================================================= */ 
/* ***********************************************************
   init
/* ***********************************************************
   reads in data, initial par values */
void  init()
{
  int i, j, code;
  double r,R,x;
  FILE *infile;

/* -----------------------------------------------------------
   read in sample size and dim
/* ----------------------------------------------------------- */
  openIn("init.mdp");
  scanInt("n", &mdp.n_obs);
  messinttab("nobs", mdp.n_obs);
  /* NOTE:  curently only q=1
     scanInt(" q ", &q); dimension of response */
  q = 1;
  scanInt("p", &p); /* DIM of regression (ANOVA) vector */
  /* scanInt(" r ", &r); 
     p = q*r; */
  messint("p ",p);
  mdp.p = p;
  scanInt("censoring", &mdp.censoring);


/* -----------------------------------------------------------
   allocate memory for pars & data
/* -----------------------------------------------------------  */
  mdp.Y = dvector(0,mdp.n_obs-1);
  mdp.CY = dvector(0,mdp.n_obs-1); /* censoring time */
  mdp.delta = ivector(0,mdp.n_obs-1); /* censoring indicator */
  mdp.d = dmatrix(0,mdp.n_obs-1,0,p-1);
  mdp.V = dvector(0,mdp.n_obs-1);
  mdp.Vj = darray_3(0,mdp.n_obs-1);
  for(j=0;j<mdp.n_obs;j++){
    mdp.Vj[j] = dmatrix(0,p-1,0,p-1);
  }
  mdp.mj = dmatrix(0,mdp.n_obs-1,0,p-1);
  mdp.mu = dmatrix(0,mdp.n_obs-1,0,p-1);
  mdp.mean = dvector(0,p-1);

  mdp.B = dmatrix(0,p-1,0,p-1);
  mdp.B_inv = dmatrix(0,p-1,0,p-1);
  mdp.B_cd_inv = dmatrix(0,p-1,0,p-1);
  mdp.C = dmatrix(0,p-1,0,p-1);

  mdp.aa = dvector(0,p-1);
  mdp.A_inv = dmatrix(0,p-1,0,p-1);

  /* note: count, first & upd allow for one extra class =
     "new class" */
  mdp.count = ivector(0,mdp.n_obs);
  mdp.first = ivector(0,mdp.n_obs);

  mdp.member = ivector(0,mdp.n_obs-1);
  mdp.updated = ivector(0,mdp.n_obs-1);
  mdp.prev = ivector(0,mdp.n_obs-1);
  mdp.next = ivector(0,mdp.n_obs);
  mdp.new = ivector(0,mdp.n_obs-1);

  member0 = ivector(0,mdp.n_obs-1);

  m_aux = dvector(0,mdp.n_obs-1);
  s_aux = dvector(0,mdp.n_obs-1);

/* -----------------------------------------------------------
   read in hyperpars and initial values
/* ----------------------------------------------------------- */

  scanInt("niter",&n_iter);
  scanInt("m_prior", &m_prior);
  scanInt("B_prior", &B_prior);
  scanInt("S_prior", &S_prior);
  scanInt("alpha_prior", &alpha_prior);
  scanInt("ndiscard",&n_discard);
  scanInt("nreinit",&n_reinitialize);
  scanInt("npi",&n_pi);
  scanInt("npredupdate",&n_predupdate);
  scanInt("nprintallpars",&n_printallpars);
  scanInt("ny0", &ny0_print);
  scanInt("verbose", &mdp.verbose);
  scanInt("px", &mdp.px);
  scanInt("s", &mdp.s);

  /* read and invert S */
  scanDouble("S-init", &mdp.S);

  scanInt("q",&mdp.q);
  scanDouble("R", &mdp.R);

  /* read and invert B */
  scanDoubleMatrix("B-init",mdp.B,p,p);
  d_invert(p,mdp.B,mdp.B_inv);
  d_chol_decomp(p,mdp.B_inv,mdp.B_cd_inv);

  scanInt("c",&mdp.c);
  scanDoubleMatrix("C",mdp.C, p,p);
  scanDoubleArray("m-init",mdp.mean, p);
  scanDoubleArray("aa", mdp.aa, p);
  scanDoubleMatrix("A-inv", mdp.A_inv, p,p);
  scanDouble("alpha-init",&mdp.alpha);
  scanDouble("a0", &mdp.a0);
  scanDouble("b0", &mdp.b0);
  if (B_prior == 1) mdp.tau = 1.0;


/* -----------------------------------------------------------
   initial config
/* ----------------------------------------------------------- */
  scanInt("k0", &n_class0);
  if (n_class0 == 1)
    for(i=0;i<mdp.n_obs;i++)
      member0[i] = 0;
  else if (n_class0 == mdp.n_obs)
    for(i=0;i<mdp.n_obs;i++)
      member0[i] = i;
  else if (n_class0 == -2){
    for(i=0;i<mdp.n_obs;i++)
      member0[i] = i%2;
    n_class0 = 2;
  }
  else
    scanIntArray("member0", member0, mdp.n_obs);

  closeIn();

/* -----------------------------------------------------------
   read in data
/* ----------------------------------------------------------- */

  infile = openIn("data.mdp");
  if (mdp.censoring==0){/* no censoring */
    scanDoubleArray(0,mdp.Y, mdp.n_obs);
    i_assgn_r(mdp.delta,1,mdp.n_obs); /* mark all as observed */
  } else {
    for(i=0;i<mdp.n_obs;i++){
      code = fscanf(infile," %lf %lf %lf", &mdp.Y[i], &x, &mdp.CY[i]);
      mdp.delta[i]=x; /* careful, delta is integer ... */
      if ( (mdp.delta[i]==0) & (mdp.Y[i]<mdp.CY[i]) ){
	printf("\n*** Error: censored (delta[%d]=0) observation has", i);
	printf("\n    (Y[%d]=%4.2f) < (C[%d]=%4.2f) ??.\n", i,mdp.Y[i],i,mdp.CY[i]);
	exit(-1);
      }
    }
  }
  closeIn();
  openIn("data-d.mdp");
  scanDoubleMatrix(0,mdp.d, mdp.n_obs, p);

  /* set r.v. seeds */
  setseed(is1,is2);

  closeIn();
}
/* ***********************************************************
   open output files
/* *********************************************************** */
void init_output()
{
  parFile = openOut("par.mdp");
  muFile = openOut("mu.mdp");
  VFile = openOut("V.mdp");
  SFile = openOut("S.mdp");
  BFile = openOut("B.mdp");
  meanFile = openOut("mean.mdp");
  fclose(muFile); 
  fclose(SFile); 
  fclose(meanFile); 
  fclose(VFile);
  fclose(BFile);
  fclose(parFile);
}

/* ***********************************************************
   initial configuration
/* *********************************************************** */
void init_config()
{
  int i, k, ct;

  /* create n empty groups */
  mdp.n_class = 0;
  for (k=0;k<mdp.n_obs;k++){
    mdp.count[k] = 0;
    mdp.first[k] = -1;
    mdp.updated[k] = 0;
    mdp.new[k] = 0;
    sample_V(k);
  }

  /* put obs's in groups */
  for (i = 0; i < mdp.n_obs; i++){
    k = member0[i];
    add_class(i,k);
  }

  /* sample mu's */
  for(k=0;k<mdp.n_obs;k++){
    sample_mu(k);
  }
  check_class();
}

/* *************************************************
   reinitialize_config
 ************************************************* */
void reinitialize_config()
{
  int i;

  /* remove old config */
  for (i=0;i<mdp.n_obs;i++)
    take_out(i);
  init_config();
}

/* ***********************************************************
   Gibbs Sampler
/* *********************************************************** */
void gibbs()
{
  int time, last_init, k;

  last_init = 0;
  for (iter = 0; iter < n_iter; iter++){
    if ( (iter % n_reinitialize == 0) & (iter > 0)){
      reinitialize_config();
      last_init = iter;
    }
    if (mdp.verbose > -1) print_pars(iter);
    if (S_prior) sample_S();
    mdp.kmax = mdp.n_class;
    sample_config();
    for(k=0; k<mdp.n_class; k++){
      sample_mu(k);
      sample_V(k);
    } 
    if ( (alpha_prior==1) | 
	 ((alpha_prior==2) & (iter > 0.5*n_discard)) )
      sample_alpha(); 
    if (B_prior) sample_B();
    if (m_prior) sample_m();
    if (iter-last_init > n_discard){
      update_nclass();
      if (iter % n_pi == 0){
	print_pi();
	write_pars(iter);
	write_nclass();
      }
      if (iter % n_predupdate == 0) /* update predictive */
	pred_update(iter, mdp);
    }
  }
}

/* ***********************************************************
   finish
/* *********************************************************** */
void finish()
{
     print_pi();
     pred_finish(mdp);
}		   
  
/* ============================================================= 
   Conditional draws
/* ============================================================= */ 
	  
/* ***********************************************************
   sample config
/* ***********************************************************
   does the Gibbs over configurations
   Pr(y[i] in class j | ...) proportional to
   for j=0,..n_class-1:
      count[j]*MV-N(y[i]; mut[j], Vt[j])
      mut[j] and Vt[j] computed in make_mut_Vt
   for j=n_class (i.e. new class): same expression with 
                                   count[n_class] = alpha
*/
void sample_config()
{
  int l, k, k_new, ind, s1;
  double u,*qv,      /* prob of pi[ind] = pi*[k] */
    logpr,mx,m,s; /* mean and var of marginal */
  int si, j;
  FILE *pltfile;

  /* allocate memory */
  qv = dvector(0,mdp.n_obs);

  for (ind = 0; ind < mdp.n_obs; ind++){
    u = runif();
    if (mdp.count[mdp.member[ind]] == 1)
      /* toss coin with prob 1/k */
      if (u > 1.0/(1.0*mdp.n_class)) continue;
    if (mdp.member[ind]==-1){
      printf("\n*** clusterless point %d! \n",ind);
      exit(1);
    }
    take_out(ind);

    /* compute prob pi[ind] = pi*[k]       */
    for(l=0;l<mdp.n_class;l++){
      if (mdp.updated[l]==0)
	setup_mu(l,0);
      logpr = loglik(ind,l,&m_aux[l],&s_aux[l]);
      qv[l] = log(1.0*mdp.count[l]) + logpr; 
    }

    /* Prob of new class */
    /* qv[n_class = alpha * MV-N(y[i],m,V+tau*B)*/
    if (mdp.updated[mdp.n_class]==0)
      setup_mu(mdp.n_class,0);
    l = mdp.n_class;
    logpr = loglik(ind,l,&m_aux[l],&s_aux[l]);
    qv[l] =log(mdp.alpha/(l+1.0))+logpr;
    /* draw the index of the new class for ind */
    mx = x_max(qv,&l,mdp.n_class+1);  /* standardize probs */
    for(l=0;l <= mdp.n_class;l++)
      qv[l]= exp(qv[l]-mx);
    multinomial(1,mdp.n_class+1,qv,&k_new);     
    if (mdp.delta[ind]==0) /* censored */
      mdp.Y[ind] = norm_trunc_left(m_aux[k_new],s_aux[k_new],mdp.CY[ind]);
    add_class(ind,k_new);

    if (k_new == mdp.n_class){ 
      if (k_new > mdp.kmax) mdp.kmax = k_new;
    }
  }

  /* free allocated memory */
  free_dvector(qv,0,mdp.n_obs);
  }


double loglik(int i, int l, double *m, double *s){
  /* compute likelihood for i-th observation with moments (mj,Vj) */
  double lp,F;
  int j,k,offset;

  xy(mdp.d[i],mdp.mj[l],m,p);
  *s = sqrt( xtAy(mdp.d[i],mdp.d[i],mdp.Vj[l],p) + mdp.V[l]);
  if (mdp.delta[i]==0){/* censored */
    F = cdfnorm(mdp.CY[i],*m,*s);
    lp = log(1.0-F);
  } else { /* not */
    lp = d_normal_dens(mdp.Y[i],*m,*s); 
  }
  return(lp);
}
	
      


/* -----------------------------------------------------------
   check class
/* ----------------------------------------------------------- */
/* checks for bogous classes */

void check_class(){
  int i,j,k;

  for(k=0;k<mdp.n_class;k++){
    for(j=0,i=mdp.first[k]; i!=-1; j++, i = mdp.next[i])
      if (j>mdp.count[k])
	printf("\n*** Error: class %d is bogous.\n", k);
  }
}

/* -----------------------------------------------------------
   take_out
/* ----------------------------------------------------------- */
void take_out(int ind)
{
  int k_old,i;
  

  /* take y[ind] out of it's current class */
  k_old = mdp.member[ind];
  mdp.count[k_old] -= 1;
  mdp.updated[k_old] = 0;
  mdp.member[ind] = -1;
  if (mdp.prev[ind] == -1)
    mdp.first[k_old] = mdp.next[ind];
  else mdp.next[mdp.prev[ind]] = mdp.next[ind];
  mdp.next[ind] = -1;
  mdp.prev[ind] = -1;

  /* if class k=k_old is empty */
  if (mdp.count[k_old] == 0){

    /* relabel class k=n_class-1 as k_old */
    mdp.updated[k_old] = mdp.updated[mdp.n_class-1];
    mdp.updated[mdp.n_class-1] = 0;
    mdp.count[k_old] = mdp.count[mdp.n_class-1];
    mdp.count[mdp.n_class-1] = 0;
    mdp.first[k_old] = mdp.first[mdp.n_class-1];
    mdp.first[mdp.n_class-1] = -1;
    mdp.new[k_old] = mdp.new[mdp.n_class-1];

    r_swap_s(&mdp.V[k_old],&mdp.V[mdp.n_class-1]);
    A_swap_B(mdp.Vj[k_old],mdp.Vj[mdp.n_class-1],mdp.p,mdp.p);
    x_swap_y(mdp.mu[k_old],mdp.mu[mdp.n_class-1],mdp.p);
    x_swap_y(mdp.mj[k_old],mdp.mj[mdp.n_class-1],mdp.p);

    /* update member-references to the relabeled class */
    for(i=mdp.first[k_old];i!=-1;i=mdp.next[i])
      mdp.member[i] = k_old;

    /* decrement n-class by one */
    mdp.n_class -= 1;
    
    check_class();

  }
}
/* -----------------------------------------------------------
   add_class
/* -----------------------------------------------------------
   adds obs[ind] to class[k]
*/
void add_class(int ind, int k)
{
  mdp.updated[k] = 0;
  mdp.count[k] += 1;
  if (mdp.count[k] == 1){
    mdp.new[k] = 1;
    mdp.n_class += 1;
  }
  mdp.member[ind] = k;
  mdp.next[ind] = mdp.first[k];
  mdp.first[k] = ind;
  if (mdp.count[k] > 1) mdp.prev[mdp.next[ind]] = ind;
  mdp.prev[ind] = -1;

  check_class();
}
/* -----------------------------------------------------------
   make_ybar
/* ----------------------------------------------------------- */
void make_yhx(int k, double *y, double *h, double **x)
{ /* make response vector, var vector and design matrix for k */
  int i,j;

  for(i=mdp.first[k],j=0;i!=-1;i=mdp.next[i],j++){
    y[j]=mdp.Y[i];
    h[j]=mdp.V[k];
    y_assgn_x(x[j],mdp.d[i],p);
  }
  if (j!=mdp.count[k]){ /* plaus check */
    printf("\n *** Error: cluster k=%d is corrupted.\n",k);
    exit(-1);
  }
}
/* -----------------------------------------------------------
   make_Sk
/* ----------------------------------------------------------- */
void make_Sk(int k, double *S)
{
  int i;
  double nk, m, z;


  nk = 1.0*mdp.count[k];
  *S = 0;
  for(i=mdp.first[k];i!=-1;i=mdp.next[i]){
    xy(mdp.d[i],mdp.mu[k],&m,p);
    z = mdp.Y[i]-m;
    *S = *S+z*z;
  }
}


/* *************************************************
   setup_mu
 ************************************************* */
void setup_mu(int k, int draw)
{
  double **Vj_cd, **Vj_inv_cd, **x, *y, *h, nk;
  double **Vj;


  /* allocate memory for the auxilary arrays */
  Vj_cd = dmatrix(0,p-1,0,p-1);
  Vj_inv_cd = dmatrix(0,p-1,0,p-1);
  y = dvector(0,mdp.count[k]);
  h = dvector(0,mdp.count[k]);
  x = dmatrix(0,mdp.count[k],0,p);

  Vj = mdp.Vj[k];

  if (mdp.count[k] == 0){
    B_assgn_A(Vj,mdp.B,mdp.p);
    y_assgn_x(mdp.mj[k],mdp.mean,p);
  }
  else{
    /* get pars for conditional posterior p(mu|...) */
    make_yhx(k,y,h,x);
    bayes_regr(mdp.mj[k], Vj_cd,Vj_inv_cd,
	       y,h,x,
	       mdp.mean, mdp.B_inv, mdp.count[k],p,1);
    ABt(Vj_cd,p,p,Vj_cd,p,p,Vj);
  }


  /* if draw, then sample mu[k] */
  if (draw==1)
    mvnS_rand(p,mdp.mj[k],Vj,mdp.mu[k]);

  mdp.updated[k] = 1;

  /* free memory */
  free_dmatrix(Vj_cd,0,p-1,0,p-1);
  free_dmatrix(Vj_inv_cd,0,p-1,0,p-1);
  free_dvector(y,0,mdp.count[k]);
  free_dvector(h,0,mdp.count[k]);
  free_dmatrix(x,0,mdp.count[k],0,p);
}

/* ***********************************************************
   sample_mu
/* ***********************************************************
   sample mu[k] from N(mu; m(y), T)
   where {m(y),T} = T*(nj*V-inv[k]*ybar[k] + B-inv/tau*m),
                    T-inv = nj*V-inv[k] + B-inv/tau.
*/
void sample_mu(int k)
{

  setup_mu(k,1);
}

/* ***********************************************************
   sample_V
   *********************************************************** */
/* V[k] from G(1/V[k]; (s+n[k])/2, 1/2*(sS+S[k])),
   where S[k] = sum (y[i]-x[i]'mu[k])^2
*/
void sample_V(int k)
{ double Sk, Vinv, s1, S1;

  /* Sk = sum (y[i]-mu[j])(y[i]-mu[j]) over member[i]=j */
  make_Sk(k,&Sk);

  /* S1 = s*S+Sk */
  S1 = 0.5*(mdp.s*mdp.S+Sk);
  s1 = 0.5*(mdp.s+mdp.count[k]);
  Vinv = gamdev_ab(s1,S1);
  mdp.V[k] = 1.0/Vinv;
}
    


/* ***********************************************************
   sample S
  *********************************************************** 
   prior:  S ~ Ga(q/2, 1/2*q/R), i.e. E(S)=R
   likl:   1/V[k] ~ Ga(s/2, 1/2*sS)
   post:   S|V    ~ Ga[1/2*(q+s*nk), 1/2(q/R + \sum s/V[k])]
*/
void sample_S()
{
  int k;

  double s,r1,R1,Vinv,S1;
  
  s = 1.0*mdp.s;
  r1 = 0.5*(mdp.q+s*mdp.n_class);
  R1 = 0.5*(mdp.q/mdp.R);
  /* loop over all k classes, add up Vbar_inv and mu */
  for(k=0;k<mdp.n_class;k++){
    R1 += 0.5*s/mdp.V[k];
  }
  /* generate S from Wish(q+s*k, U) */
  mdp.S = gamdev_ab(r1,R1);
  return;
}

/* ***********************************************************
   sample m
/* *********************************************************** 
/* m from p(m|..) = N(m; m(mubar), T)
   m(ybar[k]) = T*(B-inv/tau*k*mubar + A-inv*a)
   T-inv = B-inv/tau*k + A-inv
*/
void sample_m()
{
  double *mubar, K;
  int i;

  K = mdp.n_class; /* need double */

  /* allocate memory */
  mubar = dvector(0,p-1);
  
  /* make mubar = 1/K sum mu[k] */
  x_zero(mubar,p); /* initialize mubar=0 */
  for(i=0;i<mdp.n_class;i++){
    x_plus_y(mubar,mdp.mu[i],mubar,p);
  }
  x_div_r(mubar,K,mubar,p);

  /* compute and sample posterior on m */
  nn_bayes_rand(1.0,mdp.A_inv, mdp.aa,
		mdp.tau/K, mdp.B_inv, mubar,
		mdp.mean, p);

  /* free memory */
  free_dvector(mubar,0,p-1);
}

/* ***********************************************************
   draw B
/* *********************************************************** */
/* dummy */
/* sample B from p(B|..) =
   W(1/B; k+c, [cC + 1/tau*sum (mu[k]-m)(mu[k]-m)']^-1)
*/
void sample_B()
{ double **Sk, **S1, **S1_inv_cd,
    **B_cd, *z;
  int k;

  /* allocate memory */
  B_cd = dmatrix(0,p-1,0,p-1);
  Sk = dmatrix(0,p-1,0,p-1);
  S1 = dmatrix(0,p-1,0,p-1);
  S1_inv_cd = dmatrix(0,p-1,0,p-1);
  z = dvector(0,p-1);
  
  /* Sk = sum (mu[k]-m)(mu[k]-m)' */
  A_zero(Sk,p);
  for(k=0;k<mdp.n_class;k++){
    x_min_y(mdp.mu[k],mdp.mean,z,p);
    A_plus_rxxt(Sk,1.0,z,Sk,p);
  }

  /* S1 = s*S+Sk */
  rA_plus_sB(mdp.c*1.0, mdp.C, 1.0, Sk, S1, p);
  cdinv(p,S1,S1_inv_cd);
  wishart(mdp.c+mdp.n_class, p, S1_inv_cd, mdp.B_inv);
  
  /* compute B, B_cd, B_inv_cd etc. */
  d_chol_decomp(p,mdp.B_inv, mdp.B_cd_inv);
  d_inv_triang(p,mdp.B_cd_inv,B_cd);
  ABt(B_cd,p,p,B_cd,p,p,mdp.B);

  /* error check */
  if (isnan(mdp.B[0][0])){
    error("sample_B", "B is NaN.", 1);
  }

  /* allocate memory for the auxilary arrays */
  free_dmatrix(B_cd,0,p-1,0,p-1);
  free_dmatrix(S1,0,p-1,0,p-1);
  free_dmatrix(Sk,0,p-1,0,p-1);
  free_dmatrix(S1_inv_cd,0,p-1,0,p-1);
  free_dvector(z,0,p-1);
}



/* ***********************************************************
   sample alpha
/* *********************************************************** 
   draw alpha and eta from:
   alpha from     pi*Gamma(a0+n_class, b0-log(eta)) +
              (1-pi)*Gamma(a0+n_class-1, b0-log(eta)),
         where pi = (a+n-class-1)/(a0+n_class-1 + n*b0 - n*log(eta)).
   eta   from Beta(alpha+1,n).
*/
void sample_alpha()
{
     double a,  b, pi;       /* parameters in Gamma call */
     double u;

      dbg = 0;
      if (dbg) message("sample alpha");

     /* draw eta */
     mdp.eta = betadev(mdp.alpha+1.0, mdp.n_obs*1.0);

     /* compute pars for gammas for alpha */
     b = mdp.b0 - log(mdp.eta);
     pi = (mdp.a0 + mdp.n_class - 1)/
       (mdp.a0 + mdp.n_class -1 + mdp.n_obs*(mdp.b0 - log(mdp.eta)));
     duniform(1,&u);
     a = (u<pi) ? mdp.a0+mdp.n_class : mdp.a0+mdp.n_class-1;

     /* draw alpha with prob pi from Gamma(a1,b) else Gamma(a2,b) */
     mdp.alpha = gamdev_ab(a,b);

     if (dbg){
       messdouble("alpha:",mdp.alpha);
     }
     dbg=0;

}

/* ============================================================= 
   printing routines
/* ============================================================= */ 

/* ***********************************************************
   print pi
/* ***********************************************************
   prints current pi's */
void print_pi()
{
  int k, i,j;
     
     muFile = openAppend("mu.mdp");
     for (k = 0; k < mdp.n_class; k++){
       fprintf(muFile," %4d",iter);
       fprintf(muFile,"\t %4d \t", mdp.count[k]);
       fwriteDoubleArray(muFile, mdp.mu[k], 1, p);
     }
     fprintf(muFile,"\n");
     fclose(muFile);

     VFile = openAppend("V.mdp");
     for (k = 0; k < mdp.n_class; k++){
       fprintf(VFile,"\n %4d",iter);
       fprintf(VFile,"\t %4d ", mdp.count[k]);
       fprintf(VFile,"%5.3f\n ",mdp.V[k]);
     }
     fclose(VFile);

     if (mdp.censoring==1){
       y0file=openAppend("y0.mdp");
       for(i=0,j=0; (i<mdp.n_obs) & (j< ny0_print); i++){
	 if(mdp.delta[i]==1) /* event */	
	   continue;
	 j++;
	 fprintf(y0file," %5.3f ", mdp.Y[i]);
       }
       fprintf(y0file," \n");
       fclose(y0file);
     }
       
}

/* ***********************************************************
   print pars
/* ***********************************************************
   print current tau,alpha,eta */
void print_pars(int time)
{
  static int linect=20;
  int i,k, ct[3], k1[3];
  
  if (mdp.verbose < 2){
    printf("%d:%d ",iter,mdp.n_class);
    if (iter % 10 == 0) printf("\n");
    return;
  }

  /* print header line if linect == 20 */
  if (linect == 20){
    if (B_prior) 
      printf("\n\n\n %5s %16s %5s ", "iter","k (count1,2,3)  ","alpha");
    else
      printf("\n\n\n %5s %16s %5s %5s %5s ", "iter","k (count1,2,3)  ",
	     "alpha", "tau", "S");
    for(i=0;i<p;i++) printf("mean[%1d](B)    ",i);
    if (mdp.verbose){
      for(k=0;(k<3)&(k<mdp.n_class);k++) 
	for(i=0;i<p;i++) printf("mu[%1d][%1d]    ",k,i);
      printf
	("\n----------------------------------------------------");
      printf(
	 "-------------------------- \n");
      linect = 0;
    }
    printf("\n");
  }
  /* find largest 3 classes */
  k1[0] = k1[1] = k1[2] = 0;
  ct[0] = ct[1] = ct[2] = 0;
  for(k=0; k<mdp.n_class; k++){
    if (mdp.count[k] > ct[0]){
      k1[2] = k1[1];
      ct[2] = ct[1];
      k1[1] = k1[1];
      ct[1] = ct[0];
      k1[0] = k;
      ct[0] = mdp.count[k];
    }
    else if (mdp.count[k] > ct[1]){
      k1[2] = k1[1];
      ct[2] = ct[1];
      k1[1] = k;
      ct[1] = mdp.count[k];
    }
    else if (mdp.count[k] > ct[2]){
      k1[2] = k;
      ct[2] = mdp.count[k];
    }
  }

  /* print line */
  if (B_prior)
    printf(  "%5d %3d(%3d %3d %3d) %5.2f ",
	   iter,mdp.n_class,ct[0], ct[1], ct[2],
	   mdp.alpha);
  else
    printf(  "%5d %3d(%3d %3d %3d) %5.2f %5.2f %5.2f",
	   iter,mdp.n_class,ct[0], ct[1], ct[2],
	     mdp.alpha, mdp.tau, mdp.S);
  for(i=0;i<p;i++) printf("%4.1f(%4.1f) ", mdp.mean[i], 
			   sqrt(mdp.B[i][i]));
  if (mdp.verbose){
    for(k=0;(k<3)&(k<mdp.n_class);k++){
      if (mdp.new[k1[k]] == (iter-1)) printf("* ");
      else printf("| ");
      for(i=0;i<p;i++) printf("%4.1f ",
			      mdp.mu[k1[k]][i]);
      printf(" (%4.1f)",sqrt(mdp.V[k1[k]]));
    }
  }
  printf("\n");
  linect += 1;
}


int print_allpars()
{
  int k;

  messdoublematrix2("B ",mdp.B,p,p);
  messdouble("S ",mdp.S);
  messdoublevec("m ",mdp.mean,p);
  for(k=0;k<mdp.n_class;k++){
    messdoublevec("mu[i] ",mdp.mu[k],p);
    messdouble("V[j] ",mdp.V[k]);
  }
  printf("Note: sampling seq is \t  S,config,mu,V,B,m.\n");
  printf("------------------------------------------------------- \n");
}

/* ***********************************************************
   write pars
/* ***********************************************************
   print current tau,alpha,eta */
void write_pars(int time)
{
  parFile = openAppend("par.mdp");
  fprintf(parFile,"%4d  %6.3f %6.3f %4d %8.3e\n",
	  iter, mdp.alpha,mdp.eta, mdp.n_class, mdp.S);
  fclose(parFile);
  meanFile = openAppend("mean.mdp");
  fprintf(meanFile,"%4d \t", iter);
  fwriteDoubleArray(meanFile, mdp.mean, 1, p);
  fclose(meanFile);
  if (B_prior){
    BFile = openAppend("B.mdp");
    fprintf(BFile,"%4d \t", iter);
    fwriteDoubleMatrix2(BFile, mdp.B, p,p);
    fclose(BFile);
  }
}


/* ============================================================= 
  ergodic avges for group numbers
/* ============================================================= */ 

int nclass_initialized=0; /* indicator whether nclass is initialized */
double *nclass; /* accumulate number of clusters */
int max_nclass = 0;

/* ***********************************************************
   update_nclass
/* *********************************************************** */
int update_nclass()
{
  int j;

  if (nclass_initialized == 0){
    nclass = dvector(0,mdp.n_obs);
    for (j=0;j<=mdp.n_obs;j++)
      nclass[j] = 0.0;
    nclass_initialized = 1;
  }

  max_nclass = (mdp.n_class > max_nclass) ? mdp.n_class : max_nclass;
  nclass[mdp.n_class] += 1.0;
}


int write_nclass()
{
  int j;
  FILE *pjfile;

  pjfile = openOut("pj.mdp");
  fwriteDoubleArray(pjfile, nclass, 10, max_nclass);
  fclose(pjfile);
}
