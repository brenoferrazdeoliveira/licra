#include "../licra.h"

void ic(double *Q_00, double *Q_11, double *Q_01, double *Q_02, double *Q_12){
	int i;
	double n[3];

	const gsl_rng_type*W;
	gsl_rng*w;
	gsl_rng_env_setup();
	W= gsl_rng_default;
	w= gsl_rng_alloc(W);

	for(i= 0; i< (Nx*Ny*Nz); i++){
		gsl_ran_dir_3d(w, &n[0], &n[1], &n[2]);
		Q_00[i]= (1.0/tilde)*(0.5*S_eq*(3.0*n[0]*n[0]-1.0));
		Q_11[i]= (1.0/tilde)*(0.5*S_eq*(3.0*n[1]*n[1]-1.0));
		Q_01[i]= (1.0/tilde)*(0.5*S_eq*(3.0*n[0]*n[1]));
		Q_02[i]= (1.0/tilde)*(0.5*S_eq*(3.0*n[0]*n[2]));
		Q_12[i]= (1.0/tilde)*(0.5*S_eq*(3.0*n[1]*n[2]));
	}
	gsl_rng_free(w);
}
