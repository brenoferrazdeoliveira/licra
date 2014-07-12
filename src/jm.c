#include "../licra.h"

#define ERROr1 printf("cannot open file dat/bt-%d.dat \n", t); exit(1);
#define ERROr2 printf("cannot open file dat/jm-%d.dat \n", t); exit(1);

int main(int argc, char **argv){

	int x, y, z, i, j, k, l, o, p, N, t;
	double Ne, dol, If;	
	/* dol -> d over lambda_0 */
	double *beta, *theta;	
	double RL[2][2], RR[2][2];
	double Px[2][2]= {{1, 0}, {0, 0}, };
	double Py[2][2]= {{0, 0}, {0, 1}, };
	double complex Bi[2][2], J_all[2][2], E_f[2][2];
	double complex J[Nz][2][2];
	FILE *jm, *bt;
	char name_jm[64], name_bt[64];

	if((beta=  (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}
	if((theta= (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}

	for(t= 0; t<= NF; t++){
		sprintf(name_bt, "dat/bt-%d.dat", t);
		if((bt= fopen(name_bt, "r"))==NULL){ERROr1};
		for(N= 0; N< NxNyNz; N++){
			fscanf(bt, "%lf %lf\n", &beta[N], &theta[N]);
		}
		fclose(bt);

		sprintf(name_jm, "dat/jm-%d.dat", t);
		if((jm= fopen(name_jm, "w"))==NULL){ERROr2};
		for(y= 0; y< Ny; y++){
			for(x= 0; x< Nx; x++){
				If= 0.0;
				for(dol= 0.0; dol< 1.0; dol+=0.05){
					for(z= 0; z< Nz; z++){
						N= (z*Ny+y)*Nx+x;
			
						RL[0][0]=  cos(beta[N]);
						RL[0][1]= -sin(beta[N]);
						RL[1][0]=  sin(beta[N]);
						RL[1][1]=  cos(beta[N]);
			
						Ne= ne*no/(sqrt(ne*ne*cos(theta[N])*cos(theta[N])+\
							  no*no*sin(theta[N])*sin(theta[N])));
						Bi[0][0]= cos(2.0*M_PI*dol*Ne) - sin(2.0*M_PI*dol*Ne)*I;
						Bi[0][1]= 0.0+0.0*I;
						Bi[1][0]= 0.0+0.0*I;
						Bi[1][1]= cos(2.0*M_PI*dol*no) - sin(2.0*M_PI*dol*no)*I;
					
						RR[0][0]=  cos(beta[N]);
						RR[0][1]=  sin(beta[N]);
						RR[1][0]= -sin(beta[N]);
						RR[1][1]=  cos(beta[N]);
					
						for(i= 0; i< 2; i++){
							for(j= 0; j< 2; j++){
								J[z][i][j]= 0.0 + 0.0*I;
								for(k= 0; k< 2; k++){
									for(l= 0; l< 2; l++){
										J[z][i][j]+= RL[i][k]*Bi[k][l]*RR[l][j];
									}
								}
							}
						}
					}
#if Nz== 1
					for(o= 0; o< 2; o++){
						for(p= 0; p< 2; p++){
							J_all[o][p]= J[0][o][p];
						}
					}
#elif Nz> 1
					for(z= 0; z< (Nz-1); z++){
						for(i= 0; i< 2; i++){
							for(j= 0; j< 2; j++){
								J_all[i][j]= 0.0 + 0,0*I;
								for(k= 0; k< 2; k++){
									J_all[i][j]+= J[z][i][k]*J[z+1][k][j];
								}
								for(o= 0; o< 2; o++){
									for(p= 0; p< 2; p++){
										J[z][i][j]= J_all[i][j];
									}
								}
							}
						}
					}
#endif
					for(i= 0; i< 2; i++){
						for(j= 0; j< 2; j++){
							E_f[i][j]= 0.0 + 0.0*I;
							for(k= 0; k< 2; k++){
								for(l= 0; l< 2; l++){
									E_f[i][j]+= Px[i][k]*J_all[k][l]*Py[l][j];
								}
							}
						}
					}
					If+= cabs(E_f[0][0]+E_f[0][1])*cabs(E_f[0][0]+E_f[0][1])+cabs(E_f[1][0]+E_f[1][1])*cabs(E_f[1][0]+E_f[1][1]);
				}
				fprintf(jm, "%.2e ", If);
		}
			fprintf(jm, "\n");
		}
		fclose(jm);
	}
	free(beta);
	free(theta);
}
