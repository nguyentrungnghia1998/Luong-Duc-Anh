#include "lbm_D2Q9.h"


void D2Q9_slip_boundary_condition(double ***f, double ***ft)
{
	// ym-plane
	#pragma omp parallel for
	for(int i=1; i<im+1; i++){
		f[i][1][1]=ft[i][1][3];
		f[i][1][4]=ft[i][1][5];
		f[i][1][6]=ft[i][1][8];
	}
	// yp-plane
	#pragma omp parallel for
	for(int i=1; i<im+1; i++){
		f[i][jm][3]=ft[i][jm][1];
		f[i][jm][5]=ft[i][jm][4];
		f[i][jm][8]=ft[i][jm][6];
	}
}
