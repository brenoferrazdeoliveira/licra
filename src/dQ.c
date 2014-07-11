#include "../licra.h"

void dQ(int i, int j, int k, int N, double trQ2,
        double *trlessQ2, double *dQdt, 
        double *Q_00, double *Q_11, double *Q_01, double *Q_02, double *Q_12){
	int ip1, im1, jp1, jm1, kp1, km1;
	double D2x[5], D2y[5], D2z[5];
	double DxDy[5], DxDz[5], DyDz[5], H[5];

	ip1= (i+1)%Nx;
	im1= (i-1+Nx)%Nx;
	D2x[0]= (Q_00[(k*Ny+j)*Nx+ip1]-2.0*Q_00[N]+Q_00[(k*Ny+j)*Nx+im1]);
	D2x[1]= (Q_11[(k*Ny+j)*Nx+ip1]-2.0*Q_11[N]+Q_11[(k*Ny+j)*Nx+im1]);
	D2x[2]= (Q_01[(k*Ny+j)*Nx+ip1]-2.0*Q_01[N]+Q_01[(k*Ny+j)*Nx+im1]);
	D2x[3]= (Q_02[(k*Ny+j)*Nx+ip1]-2.0*Q_02[N]+Q_02[(k*Ny+j)*Nx+im1]);
	D2x[4]= (Q_12[(k*Ny+j)*Nx+ip1]-2.0*Q_12[N]+Q_12[(k*Ny+j)*Nx+im1]);
	
	jp1= (j+1)%Ny;
	jm1= (j-1+Ny)%Ny;
	D2y[0]= (Q_00[(k*Ny+jp1)*Nx+i]-2.0*Q_00[N]+Q_00[(k*Ny+jm1)*Nx+i]);
	D2y[1]= (Q_11[(k*Ny+jp1)*Nx+i]-2.0*Q_11[N]+Q_11[(k*Ny+jm1)*Nx+i]);
	D2y[2]= (Q_01[(k*Ny+jp1)*Nx+i]-2.0*Q_01[N]+Q_01[(k*Ny+jm1)*Nx+i]);
	D2y[3]= (Q_02[(k*Ny+jp1)*Nx+i]-2.0*Q_02[N]+Q_02[(k*Ny+jm1)*Nx+i]);
	D2y[4]= (Q_12[(k*Ny+jp1)*Nx+i]-2.0*Q_12[N]+Q_12[(k*Ny+jm1)*Nx+i]);

	DxDy[0]= (Q_00[(k*Ny+jp1)*Nx+ip1]+Q_00[(k*Ny+jm1)*Nx+im1]-Q_00[(k*Ny+jp1)*Nx+im1]-Q_00[(k*Ny+jm1)*Nx+ip1])/4.0;
	DxDy[1]= (Q_11[(k*Ny+jp1)*Nx+ip1]+Q_11[(k*Ny+jm1)*Nx+im1]-Q_11[(k*Ny+jp1)*Nx+im1]-Q_11[(k*Ny+jm1)*Nx+ip1])/4.0;
	DxDy[2]= (Q_01[(k*Ny+jp1)*Nx+ip1]+Q_01[(k*Ny+jm1)*Nx+im1]-Q_01[(k*Ny+jp1)*Nx+im1]-Q_01[(k*Ny+jm1)*Nx+ip1])/4.0;
	DxDy[3]= (Q_02[(k*Ny+jp1)*Nx+ip1]+Q_02[(k*Ny+jm1)*Nx+im1]-Q_02[(k*Ny+jp1)*Nx+im1]-Q_02[(k*Ny+jm1)*Nx+ip1])/4.0;
	DxDy[4]= (Q_12[(k*Ny+jp1)*Nx+ip1]+Q_12[(k*Ny+jm1)*Nx+im1]-Q_12[(k*Ny+jp1)*Nx+im1]-Q_12[(k*Ny+jm1)*Nx+ip1])/4.0;

#if Nz > 1
	kp1= (k+1)%Nz;
	km1= (k-1+Nz)%Nz;
	D2z[0]= (Q_00[(kp1*Ny+j)*Nx+i]-2.0*Q_00[N]+Q_00[(km1*Ny+j)*Nx+i]);
	D2z[1]= (Q_11[(kp1*Ny+j)*Nx+i]-2.0*Q_11[N]+Q_11[(km1*Ny+j)*Nx+i]);
	D2z[2]= (Q_01[(kp1*Ny+j)*Nx+i]-2.0*Q_01[N]+Q_01[(km1*Ny+j)*Nx+i]);
	D2z[3]= (Q_02[(kp1*Ny+j)*Nx+i]-2.0*Q_02[N]+Q_02[(km1*Ny+j)*Nx+i]);
	D2z[4]= (Q_12[(kp1*Ny+j)*Nx+i]-2.0*Q_12[N]+Q_12[(km1*Ny+j)*Nx+i]);

	DxDz[0]= (Q_00[(kp1*Ny+j)*Nx+ip1]+Q_00[(km1*Ny+j)*Nx+im1]-Q_00[(kp1*Ny+j)*Nx+im1]-Q_00[(km1*Ny+j)*Nx+ip1])/4.0;
	DxDz[1]= (Q_11[(kp1*Ny+j)*Nx+ip1]+Q_11[(km1*Ny+j)*Nx+im1]-Q_11[(kp1*Ny+j)*Nx+im1]-Q_11[(km1*Ny+j)*Nx+ip1])/4.0;
	DxDz[2]= (Q_01[(kp1*Ny+j)*Nx+ip1]+Q_01[(km1*Ny+j)*Nx+im1]-Q_01[(kp1*Ny+j)*Nx+im1]-Q_01[(km1*Ny+j)*Nx+ip1])/4.0;
	DxDz[3]= (Q_02[(kp1*Ny+j)*Nx+ip1]+Q_02[(km1*Ny+j)*Nx+im1]-Q_02[(kp1*Ny+j)*Nx+im1]-Q_02[(km1*Ny+j)*Nx+ip1])/4.0;
	DxDz[4]= (Q_12[(kp1*Ny+j)*Nx+ip1]+Q_12[(km1*Ny+j)*Nx+im1]-Q_12[(kp1*Ny+j)*Nx+im1]-Q_12[(km1*Ny+j)*Nx+ip1])/4.0;

	DyDz[0]= (Q_00[(kp1*Ny+jp1)*Nx+i]+Q_00[(km1*Ny+jm1)*Nx+i]-Q_00[(kp1*Ny+jm1)*Nx+i]-Q_00[(km1*Ny+jp1)*Nx+i])/4.0;
	DyDz[1]= (Q_11[(kp1*Ny+jp1)*Nx+i]+Q_11[(km1*Ny+jm1)*Nx+i]-Q_11[(kp1*Ny+jm1)*Nx+i]-Q_11[(km1*Ny+jp1)*Nx+i])/4.0;
	DyDz[2]= (Q_01[(kp1*Ny+jp1)*Nx+i]+Q_01[(km1*Ny+jm1)*Nx+i]-Q_01[(kp1*Ny+jm1)*Nx+i]-Q_01[(km1*Ny+jp1)*Nx+i])/4.0;
	DyDz[3]= (Q_02[(kp1*Ny+jp1)*Nx+i]+Q_02[(km1*Ny+jm1)*Nx+i]-Q_02[(kp1*Ny+jm1)*Nx+i]-Q_02[(km1*Ny+jp1)*Nx+i])/4.0;
	DyDz[4]= (Q_12[(kp1*Ny+jp1)*Nx+i]+Q_12[(km1*Ny+jm1)*Nx+i]-Q_12[(kp1*Ny+jm1)*Nx+i]-Q_12[(km1*Ny+jp1)*Nx+i])/4.0;

	H[0]= D2x[0]+DxDy[2]+DxDz[3]-(1.0/3.0)*(D2x[0]+D2y[1]-D2z[0]-D2z[1]+2.0*DxDy[2]+2.0*DxDz[3]+2.0*DyDz[4]);
	H[1]= DxDy[2]+D2y[1]+DyDz[4]-(1.0/3.0)*(D2x[0]+D2y[1]-D2z[0]-D2z[1]+2.0*DxDy[2]+2.0*DxDz[3]+2.0*DyDz[4]);
	H[2]= 0.5*(D2x[2]+DxDy[1]+DxDz[4]+DxDy[0]+D2y[2]+DyDz[3]);
	H[3]= 0.5*(D2x[3]+DxDy[4]-DxDz[0]-DxDz[1]+DxDz[0]+DyDz[2]+D2z[3]);
	H[4]= 0.5*(DxDy[3]+D2y[4]-DyDz[0]-DyDz[1]+DxDz[2]+DyDz[1]+D2z[4]);
	
	dQdt[0]= -Lambda*((tau+2.0*trQ2)*Q_00[N]-3.0*trlessQ2[0]-L1_tilde*(D2x[0]+D2y[0]+D2z[0])-L2_tilde*H[0]);
	dQdt[1]= -Lambda*((tau+2.0*trQ2)*Q_11[N]-3.0*trlessQ2[1]-L1_tilde*(D2x[1]+D2y[1]+D2z[1])-L2_tilde*H[1]);
	dQdt[2]= -Lambda*((tau+2.0*trQ2)*Q_01[N]-3.0*trlessQ2[2]-L1_tilde*(D2x[2]+D2y[2]+D2z[2])-L2_tilde*H[2]);
	dQdt[3]= -Lambda*((tau+2.0*trQ2)*Q_02[N]-3.0*trlessQ2[3]-L1_tilde*(D2x[3]+D2y[3]+D2z[3])-L2_tilde*H[3]);
	dQdt[4]= -Lambda*((tau+2.0*trQ2)*Q_12[N]-3.0*trlessQ2[4]-L1_tilde*(D2x[4]+D2y[4]+D2z[4])-L2_tilde*H[4]);
#elif Nz == 1
	H[0]= D2x[0]+DxDy[2]-(1.0/3.0)*(D2x[0]+D2y[1]+2.0*DxDy[2]);
	H[1]= DxDy[2]+D2y[1]-(1.0/3.0)*(D2x[0]+D2y[1]+2.0*DxDy[2]);
	H[2]= 0.5*(D2x[2]+DxDy[1]+DxDy[0]+D2y[2]);
	H[3]= 0.5*(D2x[3]+DxDy[4]);
	H[4]= 0.5*(DxDy[3]+D2y[4]);
	
	dQdt[0]= -Lambda*((tau+2.0*trQ2)*Q_00[N]-3.0*trlessQ2[0]-L1_tilde*(D2x[0]+D2y[0])-L2_tilde*H[0]);
	dQdt[1]= -Lambda*((tau+2.0*trQ2)*Q_11[N]-3.0*trlessQ2[1]-L1_tilde*(D2x[1]+D2y[1])-L2_tilde*H[1]);
	dQdt[2]= -Lambda*((tau+2.0*trQ2)*Q_01[N]-3.0*trlessQ2[2]-L1_tilde*(D2x[2]+D2y[2])-L2_tilde*H[2]);
	dQdt[3]= -Lambda*((tau+2.0*trQ2)*Q_02[N]-3.0*trlessQ2[3]-L1_tilde*(D2x[3]+D2y[3])-L2_tilde*H[3]);
	dQdt[4]= -Lambda*((tau+2.0*trQ2)*Q_12[N]-3.0*trlessQ2[4]-L1_tilde*(D2x[4]+D2y[4])-L2_tilde*H[4]);
#endif
}
