#include "../licra.h"

#define ERROr1 printf("cannot open file dat/nl-%d.dat \n", file); exit(1);
#define ERROr2 printf("cannot open file dat/bt-%d.dat \n", file); exit(1);
#define ERROr3 printf("cannot open file dat/cp-%d.dat \n", file); exit(1);
#define ERROr4 printf("cannot open file dat/S-%d.dat  \n", file); exit(1);
#define ERROr5 printf("cannot open file dat/P-%d.dat  \n", file); exit(1);

void op(int file,
        double *Q_00, double *Q_11, double *Q_01, double *Q_02, double *Q_12){
	int i, j, N;
	double temp, d0, d1, d2, n[3], l[3], data[9];
	FILE *cp, *bt, *S, *P;
	char name_cp[64], name_bt[64], name_S[64], name_P[64];

#if Nz== 1
	N= 0;
	FILE *nl;
	char name_nl[64];
	sprintf(name_nl, "dat/nl-%d.dat", file);
	sprintf(name_bt, "dat/bt-%d.dat", file);
	sprintf(name_cp, "dat/cp-%d.dat", file);
	sprintf(name_S,  "dat/S-%d.dat",  file);
	sprintf(name_P,  "dat/P-%d.dat",  file);
	if((nl= fopen(name_nl, "w"))==NULL){ERROr1};
	if((bt= fopen(name_bt, "w"))==NULL){ERROr2};
	if((cp= fopen(name_cp, "w"))==NULL){ERROr3};
	if((S = fopen(name_S , "w"))==NULL){ERROr4};
	if((P = fopen(name_P , "w"))==NULL){ERROr5};
	for(j= 0; j< Ny; j++){
		for(i= 0; i< Nx; i++){
			data[0]= tilde*Q_00[N];
			data[1]= tilde*Q_01[N];
			data[2]= tilde*Q_02[N];
			data[3]= tilde*Q_01[N];
			data[4]= tilde*Q_11[N];
			data[5]= tilde*Q_12[N];
			data[6]= tilde*Q_02[N];
			data[7]= tilde*Q_12[N];
			data[8]= tilde*(-Q_00[N]-Q_11[N]);

			gsl_matrix_view m= gsl_matrix_view_array(data, 3, 3);
			gsl_vector *eval= gsl_vector_alloc(3);
			gsl_matrix *evec= gsl_matrix_alloc(3, 3);
			gsl_eigen_symmv_workspace * w= gsl_eigen_symmv_alloc(3);
			gsl_eigen_symmv(&m.matrix, eval, evec, w);
			gsl_eigen_symmv_free(w);
			gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

			d0= gsl_vector_get(eval, 0);
			d1= gsl_vector_get(eval, 1);
//			d2= gsl_vector_get(eval, 2);

			n[0]= gsl_matrix_get(evec, 0, 0);
			n[1]= gsl_matrix_get(evec, 1, 0);
			n[2]= gsl_matrix_get(evec, 2, 0);

			l[0]= gsl_matrix_get(evec, 0, 1);
			l[1]= gsl_matrix_get(evec, 1, 1);
//			l[2]= gsl_matrix_get(evec, 2, 1);

			gsl_vector_free(eval);
			gsl_matrix_free(evec);

			temp= atan2(n[1], n[0]);
			fprintf(nl, "%d %d %e %e %e %e \n", i, j, n[0], n[1], l[0], l[1]);
			fprintf(bt, "%e %e\n", atan2(n[1], n[0]), acos(n[2]));
			fprintf(cp, "%e ", sin(2.0*temp)*sin(2.0*temp));
			fprintf(S,  "%e ", d0);
			fprintf(P,  "%e ", 2.0*d1+d0);
			N++;
		}
		fprintf(cp, "\n");
		fprintf(S,  "\n");
		fprintf(P,  "\n");
	}
	fclose(nl);
	fclose(bt);
	fclose(cp);
	fclose(S);
	fclose(P);
#elif Nz> 1
	sprintf(name_bt, "dat/bt-%d.dat", file);
	sprintf(name_S,  "dat/S-%d.dat",  file);
	sprintf(name_P,  "dat/P-%d.dat",  file);
	if((bt= fopen(name_bt, "w"))==NULL){ERROr2};
	if((S = fopen(name_S , "w"))==NULL){ERROr4};
	if((P = fopen(name_P , "w"))==NULL){ERROr5};
	for(N= 0; N< NxNyNz; N++){
		data[0]= tilde*Q_00[N];
		data[1]= tilde*Q_01[N];
		data[2]= tilde*Q_02[N];
		data[3]= tilde*Q_01[N];
		data[4]= tilde*Q_11[N];
		data[5]= tilde*Q_12[N];
		data[6]= tilde*Q_02[N];
		data[7]= tilde*Q_12[N];
		data[8]= tilde*(-Q_00[N]-Q_11[N]);

		gsl_matrix_view m= gsl_matrix_view_array(data, 3, 3);
		gsl_vector *eval= gsl_vector_alloc(3);
		gsl_matrix *evec= gsl_matrix_alloc(3, 3);
		gsl_eigen_symmv_workspace * w= gsl_eigen_symmv_alloc(3);
		gsl_eigen_symmv(&m.matrix, eval, evec, w);
		gsl_eigen_symmv_free(w);
		gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

		d0= gsl_vector_get(eval, 0);
		d1= gsl_vector_get(eval, 1);
//		d2= gsl_vector_get(eval, 2);

		n[0]= gsl_matrix_get(evec, 0, 0);
		n[1]= gsl_matrix_get(evec, 1, 0);
		n[2]= gsl_matrix_get(evec, 2, 0);

//		l[0]= gsl_matrix_get(evec, 0, 1);
//		l[1]= gsl_matrix_get(evec, 1, 1);
//		l[2]= gsl_matrix_get(evec, 2, 1);

		gsl_vector_free(eval);
		gsl_matrix_free(evec);

		fprintf(bt, "%e %e\n", atan2(n[1], n[0]), acos(n[2]));
		fprintf(S,  "%e\n", d0);
		fprintf(P,  "%e\n", 2.0*d1+d0);
	}
	fclose(bt); 
	fclose(S);
	fclose(P);
#endif
}
