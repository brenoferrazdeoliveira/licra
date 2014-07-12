#include "../licra.h"

#define ERROr6 printf("cannot open file dat/Lxt.dat\n"); exit(1);

void sd(int m, int counter,
        double *Q_00, double *Q_11, double *Q_01, double *Q_02, double *Q_12){
	int x, y, z, N, NT, Nxh= Nx/2+1;
	double kx, ky, kz, mean_value_k2, sum_Sxk, *Sxk;
	FILE *out_L;

	if((out_L= fopen("dat/Lxt.dat", "a"))==NULL){ERROr6};

	if((Sxk= (double *)malloc(NxNyNz*sizeof(double)))==NULL){ERROr}

	fftw_complex *out_00= fftw_malloc(sizeof(fftw_complex)*(Nxh*Ny*Nz));
	fftw_complex *out_11= fftw_malloc(sizeof(fftw_complex)*(Nxh*Ny*Nz));
	fftw_complex *out_01= fftw_malloc(sizeof(fftw_complex)*(Nxh*Ny*Nz));
	fftw_complex *out_02= fftw_malloc(sizeof(fftw_complex)*(Nxh*Ny*Nz));
	fftw_complex *out_12= fftw_malloc(sizeof(fftw_complex)*(Nxh*Ny*Nz));

	fftw_plan fft_00= fftw_plan_dft_r2c_3d(Nz, Ny, Nx, Q_00, out_00, FFTW_ESTIMATE);
	fftw_plan fft_11= fftw_plan_dft_r2c_3d(Nz, Ny, Nx, Q_11, out_11, FFTW_ESTIMATE);
	fftw_plan fft_01= fftw_plan_dft_r2c_3d(Nz, Ny, Nx, Q_01, out_01, FFTW_ESTIMATE);
	fftw_plan fft_02= fftw_plan_dft_r2c_3d(Nz, Ny, Nx, Q_02, out_02, FFTW_ESTIMATE);
	fftw_plan fft_12= fftw_plan_dft_r2c_3d(Nz, Ny, Nx, Q_12, out_12, FFTW_ESTIMATE);

	fftw_execute(fft_00);
	fftw_execute(fft_11);
	fftw_execute(fft_01);
	fftw_execute(fft_02);
	fftw_execute(fft_12);

	sum_Sxk= 0.0;
	NT= 0;
	for(z= 0; z< Nz; z++){
		for(y= 0; y< Ny; y++){
			N= (z*Ny+y)*Nxh;
			Sxk[NT]= out_00[N]*conj(out_00[N])+
			         out_11[N]*conj(out_11[N])+
			         (-out_00[N]-out_11[N])*conj(-out_00[N]-out_11[N])+
			         2.0*out_01[N]*conj(out_01[N])+
			         2.0*out_02[N]*conj(out_02[N])+
			         2.0*out_12[N]*conj(out_12[N]);
			sum_Sxk+= Sxk[NT];
			NT++;
			for(x= 1; x< (Nxh-1); x++){
				N= (z*Ny+y)*Nxh+x;
				Sxk[NT]= out_00[N]*conj(out_00[N])+
				         out_11[N]*conj(out_11[N])+
				         (-out_00[N]-out_11[N])*conj(-out_00[N]-out_11[N])+
				         2.0*out_01[N]*conj(out_01[N])+
				         2.0*out_02[N]*conj(out_02[N])+
				         2.0*out_12[N]*conj(out_12[N]);
				sum_Sxk+= Sxk[NT];
				NT++;
			}
			N= (z*Ny+y)*Nxh+Nxh-1;
			Sxk[NT]= out_00[N]*conj(out_00[N])+
			         out_11[N]*conj(out_11[N])+
			         (-out_00[N]-out_11[N])*conj(-out_00[N]-out_11[N])+
			         2.0*out_01[N]*conj(out_01[N])+
			         2.0*out_02[N]*conj(out_02[N])+
			         2.0*out_12[N]*conj(out_12[N]);
			sum_Sxk+= Sxk[NT];
			NT++;
			for(x= (Nxh-2); x> 0; x--){
				N= (z*Ny+y)*Nxh+x;
				Sxk[NT]= out_00[N]*conj(out_00[N])+
				         out_11[N]*conj(out_11[N])+
				         (-out_00[N]-out_11[N])*conj(-out_00[N]-out_11[N])+
				         2.0*out_01[N]*conj(out_01[N])+
				         2.0*out_02[N]*conj(out_02[N])+
				         2.0*out_12[N]*conj(out_12[N]);
				sum_Sxk+= Sxk[NT];
				NT++;
			}
		}
	}

	for(N= 0; N< NxNyNz; N++){
		Sxk[N]/= sum_Sxk;
	}

	mean_value_k2= 0.0;
	N= 0;
	for(z= 0; z< Nz; z++){
		if(z< (Nz/2+1)){kz= (double)z/(double)Nz;}
		else{kz= (double)(z-Nz)/(double)Nz;}
		for(y= 0; y< Ny; y++){
			if(y< (Ny/2+1)){ky= (double)y/(double)Ny;}
			else{ky= (double)(y-Ny)/(double)Ny;}
			for(x= 0; x< Nx; x++){
				if(x< Nxh){kx= (double)x/(double)Nx;}
				else{kx= (double)(x-Nx)/(double)Nx;}
				mean_value_k2+= (kz*kz+ky*ky+kx*kx)*Sxk[N];
				N++;
			}
		}
	}

#if Nz== 1
	fprintf(out_L, "%e %e\n", counter*dt, sqrt(1.0/(mean_value_k2*Nx*Ny)));
#elif Nz> 1
	fprintf(out_L, "%e %e\n", counter*dt, sqrt(1.0/(mean_value_k2*cbrt(NxNyNz)*cbrt(NxNyNz))));
#endif

	fclose(out_L);
	free(Sxk);

	fftw_destroy_plan(fft_00);
	fftw_destroy_plan(fft_11);
	fftw_destroy_plan(fft_01);
	fftw_destroy_plan(fft_02);
	fftw_destroy_plan(fft_12);

	fftw_free(out_00);
	fftw_free(out_11);
	fftw_free(out_01);
	fftw_free(out_02);
	fftw_free(out_12);

	fftw_cleanup();
}
