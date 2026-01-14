/* ===================================================================== 
   Predictives and cond expectations in MDP

   to be used with mdp.c in the same directory
/* ===================================================================== */ 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "rand.h"
#include "nrutil.h"
#include "matrix.h"
#include "vector.h"
#include "rand-nr.h" 
#include "mess.h"
#include "interface.h"
#include "ddpanova.h"
#include "predict.h"

int makegrid(double *grid, int m, double from, double to);
void write_pi(int d, double *S, int nx);

/* ***********************************************************
   data structures
/* *********************************************************** */

static FILE *zfile;
static FILE *piFile;

/* -----------------------------------------------------------
   grid for univariate pdf and exp value
/* ----------------------------------------------------------- */
struct X_GRID{
  int nx;
  int nd,              /* number rows in d */
    nd0;               /* posterior simulations for the first nd0
			  designs (of nd) are saved */
  int sim;             /* indicator for including predictve sim */
  int nsim;            /* how many each time.. */
  double **d;          /* design vector */
  int *resid;          /* indicator for d-th prediction to
			  be fit or predictive inference */
  int n_update;        /* number of times z-value has been upd */
  double *ygrid;       /* grid on x - if applicable */
  double **p,          /* grid on p values if applic */  
    **Sy,**Sy2,	       /* survival and sd. */
    **hy,**hy2;	       /* hazard and sd. */
  double **S;          /* cdf on left and right end of grid */
  int xlist;           /* indicator for evaluating 
			  pdf on mdp data */
  double **pX;         /* z values on mdp data */
};

struct X_GRID *xalloc(void)
/* allocates memory for an X_GRID */
{
  return (struct X_GRID *) malloc(sizeof(struct X_GRID));
}

                

/* ***********************************************************
   univariate predictive
/* *********************************************************** */
/* -----------------------------------------------------------
   initialize univ predicitve
/* ----------------------------------------------------------- */
struct X_GRID *px_init(struct MDP_PARS mdp)
{
  int nx, nd, nd0, xi, zi, p, sim, nsim,i;
  double lo,hi;
  struct X_GRID *px;

  p = mdp.p;

  /* read pars */
  openIn("px-init.mdp");
  scanInt("nx", &nx);
  scanInt("nd", &nd);
  scanInt("nd0", &nd0);
  scanInt("sim", &sim);
  scanInt("nsim", &nsim);
  /* alloc mem */
  px = xalloc();
  px->p = dmatrix(0,nd-1,0,nx-1);
  px->Sy = dmatrix(0,nd-1,0,nx-1);
  px->Sy2 = dmatrix(0,nd-1,0,nx-1);
  px->hy = dmatrix(0,nd-1,0,nx-1);
  px->hy2 = dmatrix(0,nd-1,0,nx-1);
  px->d = dmatrix(0,nd-1,0,p-1);
  px->ygrid = dvector(0,nx-1);
  px->resid = ivector(0,nd-1);

  /* init members */
  px->nx = nx;
  px->nd = nd;
  px->nd0 = nd0;
  px->sim = sim;
  px->nsim = nsim;

  scanDoubleMatrix("d", px->d, nd,p);
  scanIntArray("resid", px->resid, nd);
  px->n_update = 0;

  /* read x grid */
  scanDoubleArray("ygrid", px->ygrid,2); // first read in hi/lo
  lo=px->ygrid[0];
  hi=px->ygrid[1];
  makegrid(px->ygrid,nx,lo,hi);


  /* check for evaluation on mdp data */
  scanInt("xlist", &px->xlist);
  if (px->xlist==1)
    px->pX = dmatrix(0,nd-1,0,mdp.n_obs-1);
  
  closeIn();

  if(px->sim)
    zfile=openOut("z.mdp");
  if(nd0>0)
    piFile=openOut("Si.mdp");

  return px;
}

/* -----------------------------------------------------------
   print_px
/* ----------------------------------------------------------- */
/*   print E(z|xy) */
void print_px(struct X_GRID *px, struct MDP_PARS mdp)
{
  int p,i,d;
  FILE *sSfile, *shfile;
  double sSi, shi;
  p=mdp.p;

  /* print x-grid */
  openOut("px-x.mdp");
  writeDoubleArray(px->ygrid,1,px->nx);
  closeOut();

  /* print d */
  openOut("px-d.mdp");
  writeDoubleMatrix(px->d,px->nd,p);
  closeOut();

  /* write pdf on grid */
  openOut("px-p.mdp");
  writeDoubleMatrix(px->p,px->nd,px->nx);
  closeOut();

  /* write cdf on grid */
  openOut("px-Sy.mdp");
  writeDoubleMatrix(px->Sy,px->nd,px->nx);
  closeOut();
  openOut("px-hy.mdp");
  writeDoubleMatrix(px->hy,px->nd,px->nx);
  closeOut();
  /* sd's */
  sSfile= openOut("px-sSy.mdp");
  shfile= openOut("px-shy.mdp");
  for(d=0;d<px->nd;d++){
    for(i=0;i<px->nx;i++){
      sSi = sqrt(px->Sy2[d][i] - px->Sy[d][i]*px->Sy[d][i]);
      shi = sqrt(px->hy2[d][i] - px->hy[d][i]*px->hy[d][i]);
      fprintf(sSfile," %5.3e ", sSi);
      fprintf(shfile," %5.3e ", shi);
    }
  }
  fclose(sSfile);
  fclose(shfile);

  /* write pred on mdp data */
  if (px->xlist==1){
    openOut("px-pX.mdp");
    writeDoubleMatrix(px->pX,px->nd,mdp.n_obs);
    closeOut();
  }

  /* if posterior predictive draws */

}

/* -----------------------------------------------------------
   update
/* ----------------------------------------------------------- */
/*   get expected values over grid xgrid, ygrid (global??) */
     void px_update(int time, struct X_GRID *px, 
		struct MDP_PARS mdp)
{
  int i, j, k, p, d, nx, xi, n_updates;
  double y, sxx, m, sx, r, r0,r1,*z, *S, *pX, kd;
  double pi, alpha, an, *u, *zi;
  double Si, hi; /* cdf on left and right end of the grid */ 

  /* set aux pars */  
  p = mdp.p;
  nx = px->nx;
  alpha = mdp.alpha;
  kd = 1.0*mdp.n_class;
  an = mdp.alpha+1.0*mdp.n_obs;

  /* alloc mem */
  z = dvector(0,nx-1);
  S = dvector(0,nx-1);
  pX = dvector(0,mdp.n_obs);
  u = dvector(0,px->nsim);
  zi = dvector(0,p);

  /* uniform r.v. for nsim r.v. from posterior on theta */
  duniform(px->nsim, u);
  
  for(d=0;d<px->nd;d++){
    /* initialize L for V = (1+tau)*s/(s-1)*S */
    k = mdp.n_class;
    xy(px->d[d],mdp.mj[k],&m,p);  /* CHECK THAT - or use mdp.mean? */
    sxx =  xtAy(px->d[d],px->d[d],mdp.Vj[k],p);
    if (px->resid[d]==1)
      sxx = sxx + mdp.V[k];
    sx = sqrt(sxx);
    r = alpha/an;

    /* predictive if y from G0 */
    for (i = 0; i < nx; i++){
      y = px->ygrid[i];
      pi = pdfnorm(y, m, sx);
      Si = 1.0-cdfnorm(y,m,sx);
      z[i] = r*pi;
      S[i] = r*Si;
    }

    if (px->xlist==1)
      for(i=0; i<mdp.n_obs;i++){
	y = mdp.Y[i];
	pi = pdfnorm(y,m,sx);
	pX[i] = r*pi;
      }
    if ((px->sim==1) & (d==0))
      zsim(u,px->nsim,r,mdp.mj[k],mdp.Vj[k],p,zi);

  
    /* pred if y is in one of the pxisting classes */
    for (k = 0; k < mdp.n_class; k++){
      r = 1.0*mdp.count[k]/an;
      xy(px->d[d],mdp.mj[k],&m,p);
      sxx =  xtAy(px->d[d],px->d[d],mdp.Vj[k],p);
      if (px->resid[d]==1)
	sxx = sxx + mdp.V[k];
      sx = sqrt(sxx);
      for (i = 0; i < nx; i++){
	y = px->ygrid[i];
	/* if (fabs(y-m) < 3*sx) */
	pi = pdfnorm(y, m, sx);
	z[i] += r*pi;
	Si = 1.0-cdfnorm(y,m,sx);
	S[i] += r*Si;
      }
      if (px->xlist==1)
	for(i=0; i<mdp.n_obs;i++){
	  y = mdp.Y[i];
	  pi = pdfnorm(y,m,sx);
	  pX[i] += r*pi;
	}
      if ((px->sim==1) & (d==0))
	zsim(u,px->nsim,r,mdp.mj[k],mdp.Vj[k],p,zi);
    }

    /* update px->p, Sy, hy */
    n_updates = px->n_update;	
    r0 = n_updates*1.0/(1.0+n_updates);
    r1 = 1.0-r0;
    for (i = 0; i < nx; i++){ /* pdf */
      px->p[d][i] = r0*px->p[d][i] + r1*z[i];
      Si = S[i];      /* S(y[i]) */
      if (Si <0){
	printf("\n *** Error: S[nx] does not match S[1]-sum p[i].\n");
	exit(-1);
      }
      if (Si < 0.00001) Si = 0.00001;
      hi = z[i]/Si;   /* h(y[i]) */
      px->hy[d][i] = r0*px->hy[d][i] + r1*hi;
      px->Sy[d][i] = r0*px->Sy[d][i] + r1*Si;
      px->hy2[d][i] = r0*px->hy2[d][i] + r1*hi*hi;
      px->Sy2[d][i] = r0*px->Sy2[d][i] + r1*Si*Si;
    }

    /* update px->pX */
    if (px->xlist == 1)
      for(i=0; i<mdp.n_obs;i++)
	px->pX[d][i] = (n_updates* (px->pX[d][i]) + pX[i])/(n_updates+1.0);
    if(d<px->nd0){ /* save imputed p */
      write_pi(d,S, px->nx);
    }
  }/* for d*/
  
  /* increment n_update */
  px->n_update = n_updates+1;

  /* free alloc mem */
  free_dvector(z,0,nx-1);
  free_dvector(S,0,nx-1);
  free_dvector(pX,0,mdp.n_obs);
  free_dvector(u,0,px->nsim);
  free_dvector(zi,0,p);
  return;
}        


void zsim(double *u, int nsim, double pk, 
	  double *m, double **S, int p, double *z)
{ /* generates theta[i] for each u[i] with u[i]-pk <= 0 */
  int i;
  
  x_plus_r(u,-pk,u,nsim);
  for(i=0;i<nsim;i++){
    if (u[i]<=0){ /* generate.. */
      mvnS_rand(p,m,S,z);
      fwriteDoubleArray(zfile,z,1,p);
    }
  }
}

void write_pi(int d, double *S, int nx)
{/* aux function to write out current guy.. */
  int i;
  
  fprintf(piFile,"%d ",d);
  for(i=0;i<nx;i++)
    fprintf(piFile," %5.3f ",S[i]);
  fprintf(piFile,"\n");
}

/* ===================================================================== 
   problem specific pred routines
/* ===================================================================== */ 

static struct X_GRID *px;

/* ***********************************************************
   pred_init
/* *********************************************************** */
void pred_init(struct MDP_PARS mdp)
{
  int j;

  message("pred_init");
  if (mdp.px == 1)   px = px_init(mdp);
}

/* ***********************************************************
   pred_update
/* *********************************************************** */
void pred_update(int time, struct MDP_PARS mdp)
{
  if (mdp.px == 1) px_update(time,px,mdp);
  pred_finish(mdp);
}
/* ***********************************************************
   pred_finish
/* *********************************************************** */
void pred_finish(struct MDP_PARS mdp)
{
  if (mdp.px == 1){
    print_px(px,mdp);
    if (px->sim==1) 
      fflush(zfile);
    if (px->nd0>0) 
      fflush(piFile);
  }
}

/* *************************************************
   makegrid
 ************************************************* */
int makegrid(double *grid, int m, double from, double to)
{
  int j;
  double step;

  step = (to-from)/((m-1)*1.0);
  for(j=0;j < m;j++)
    grid[j] = from+step*j;
}
