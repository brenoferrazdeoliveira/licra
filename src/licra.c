#include "../licra.h"

#define total (tf/dt)
#define file (int)(total/NF)

void ic(double *, double *, double *, double *, double *);
void tr(int N, double *trQ2, double *trlessQ2,
        double *, double *, double *, double *, double *);
void dQ(int i, int j, int k, int N, double trQ2, 
        double *trlessQ2, double *dQdt,
        double *, double *, double *, double *, double *);
void sd(int m, int counter,
        double *, double *, double *, double *, double *);
void op(int counter, 
        double *, double *, double *, double *, double *);

int main(int argc, char **argv){
	PRINT
	int i, j, k, l, m= 40, N, counter= 0;
	double trQ2;
	double *trlessQ2, *dQdt;
	double *Q_00,     *Q_11,     *Q_01,     *Q_02,     *Q_12;
	double *Qtemp_00, *Qtemp_11, *Qtemp_01, *Qtemp_02, *Qtemp_12;

	if((trlessQ2= (double *)malloc(5     *sizeof(double)))==NULL){ERROr}
	if((dQdt=     (double *)malloc(5     *sizeof(double)))==NULL){ERROr}
	
	if((Q_00=     (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((Q_11=     (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((Q_01=     (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((Q_02=     (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((Q_12=     (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	
	if((Qtemp_00= (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((Qtemp_11= (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((Qtemp_01= (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((Qtemp_02= (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((Qtemp_12= (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	
	ic(   Q_00, Q_11, Q_01, Q_02, Q_12);
	op(0, Q_00, Q_11, Q_01, Q_02, Q_12);
	 
	while(counter< total){
		N= 0;
		for(k= 0; k< Nz; k++){
			for(j= 0; j< Ny; j++){
				for(i= 0; i< Nx; i++){
					tr(N, &trQ2, trlessQ2,
	           Q_00, Q_11, Q_01, Q_02, Q_12);
					dQ(i, j, k, N, trQ2, trlessQ2, dQdt, 
	           Q_00, Q_11, Q_01, Q_02, Q_12);
					Qtemp_00[N]= Q_00[N]+0.5*dQdt[0];
					Qtemp_11[N]= Q_11[N]+0.5*dQdt[1];
					Qtemp_01[N]= Q_01[N]+0.5*dQdt[2];
					Qtemp_02[N]= Q_02[N]+0.5*dQdt[3];
					Qtemp_12[N]= Q_12[N]+0.5*dQdt[4];
					N++;
				}
			}
		}
		N= 0;
		for(k= 0; k< Nz; k++){
			for(j= 0; j< Ny; j++){
				for(i= 0; i< Nx; i++){
					tr(N, &trQ2, trlessQ2, 
	           Qtemp_00, Qtemp_11, Qtemp_01, Qtemp_02, Qtemp_12);
					dQ(i, j, k, N, trQ2, trlessQ2, dQdt, 
	           Qtemp_00, Qtemp_11, Qtemp_01, Qtemp_02, Qtemp_12);
					Q_00[N]= Q_00[N]+dQdt[0];
					Q_11[N]= Q_11[N]+dQdt[1];
					Q_01[N]= Q_01[N]+dQdt[2];
					Q_02[N]= Q_02[N]+dQdt[3];
					Q_12[N]= Q_12[N]+dQdt[4];
					N++;
				}
			}
		}
		counter++;
		if(counter>= (int)(exp(log(total)/100.0*m))){
			sd(m, counter, Q_00, Q_11, Q_01, Q_02, Q_12);
			m++;
		}
		if((counter%file)==0){
			op((counter)/file, Q_00, Q_11, Q_01, Q_02, Q_12);
			printf("   %.1f%%\n", (counter*100)/total);
		}
	}

	free(trlessQ2);
	free(dQdt);
	free(Q_00); free(Qtemp_00);
	free(Q_11); free(Qtemp_11);
	free(Q_01); free(Qtemp_01);
	free(Q_02); free(Qtemp_02);
	free(Q_12); free(Qtemp_12);
}
