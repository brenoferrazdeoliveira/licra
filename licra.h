#include <stdio.h>            /*************************************/
#include <stdlib.h>           /*                                   */
#include <string.h>           /*              licra                */
#include <math.h>             /*           version 3.1             */
#include <gsl/gsl_math.h>     /*                                   */
#include <gsl/gsl_rng.h>      /*           2014/Jul/07             */
#include <gsl/gsl_randist.h>  /*                                   */
#include <gsl/gsl_eigen.h>    /*          licra@fc.up.pt           */
#include <complex.h>          /*                                   */
#include <fftw3.h>            /*************************************/

#define Nx 256
#define Ny 256
#define Nz 1
#define NF 3   /* number of output files              */
#define FT 0   /* 0 -> png | 1 -> pdf - output format */

#if LC==0
	/* MBBA */
	#define PRINT printf("-> MBBA <-\n");
	#define a 0.0867            /* 10^6 J/Km^3 */
	#define B -2.12             /* 10^6 J/m^3  */
	#define C 1.74              /* 10^6 J/m^3  */
	#define L1 3.5              /* 10^-12 N    */
	#define L2 0.0              /* 10^-12 N    */
	#define dx 1.0              /* 10^-9 m     */
	#define dt 4.0              /* 10^-9 s     */
	#define tf 8e3              /* 10^-9 s     */
	#define mu_1 0.2            /* Pa s        */
	#define T (34.0-38.0)       /* (T-T*)      */
	#define ne 1.7
	#define no 1.5
#elif LC==1
	/* 5CB */
	#define PRINT printf("-> 5CB <-\n");
	#define a 0.13              /* 10^6 J/Km^3 */
	#define B -1.60             /* 10^6 J/m^3  */
	#define C 3.90              /* 10^6 J/m^3  */
	#define L1 1.6              /* 10^-12 N    */
	#define L2 4.9              /* 10^-12 N    */
	#define dx 1.0              /* 10^-9 m     */
	#define dt 2.0              /* 10^-9 s     */
	#define tf 10e3             /* 10^-9 s     */
	#define mu_1 0.3            /* Pa s        */
	#define T (26.0-34.2)       /* (T-T*)      */
	#define ne 1.7
	#define no 1.5
#endif

#define tau ((9.0*C*a*T)/(2.0*B*B))
#define S_eq (((-B+sqrt((B*B-24.0*a*T*C))))/(6.0*C))
#define L1_tilde ((9.0*C*L1)/(2.0*B*B*dx*dx))
#define L2_tilde ((9.0*C*L2)/(2.0*B*B*dx*dx))
#define tilde ((-2.0*B)/(3.0*C))
#define Lambda ((2.0*B*B*dt*1.0e-3)/(9.0*C*mu_1))

#define NxNyNz Nx*Ny*Nz
#define ERROr printf("cannot allocate memory\n"); exit(1);
