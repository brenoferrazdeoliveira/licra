#include "../licra.h"

void tr(int N, double *trQ2, double *trlessQ2, 
        double *Q_00, double *Q_11, double *Q_01, double *Q_02, double *Q_12){

	*trQ2= 2.0*(
	            Q_00[N]*Q_00[N]+ Q_11[N]*Q_11[N]+ 
	            Q_01[N]*Q_01[N]+ Q_02[N]*Q_02[N]+
	            Q_12[N]*Q_12[N]+ Q_00[N]*Q_11[N]
	           );

	trlessQ2[0]= Q_00[N]*Q_00[N]+ Q_01[N]*Q_01[N]+
	             Q_02[N]*Q_02[N]- (1.0/3.0)*(*trQ2);

	trlessQ2[1]= Q_01[N]*Q_01[N]+ Q_11[N]*Q_11[N]+
	             Q_12[N]*Q_12[N]- (1.0/3.0)*(*trQ2);

	trlessQ2[2]= Q_01[N]*Q_00[N]+ Q_11[N]*Q_01[N]+
	             Q_12[N]*Q_02[N];

	trlessQ2[3]= Q_02[N]*Q_00[N]+ Q_12[N]*Q_01[N]+
	             (-Q_00[N]-Q_11[N])*Q_02[N];

	trlessQ2[4]= Q_02[N]*Q_01[N]+ Q_12[N]*Q_11[N]+
	             (-Q_00[N]-Q_11[N])*Q_12[N];
}
