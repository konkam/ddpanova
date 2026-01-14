void pred_init(struct MDP_PARS);
void pred_update(int time, struct MDP_PARS mdp);
void pred_finish(struct MDP_PARS);
void zsim(double *u, int nsim, double pk, 
	  double *m, double **S, int p, double *z1);
