#include "lbm_D2Q9.h"


void D2Q9_domain_bc_periodic(double ***f, double ***ft)
{
	#pragma omp parallel for
	for(int j=1; j<jm+1; j++){

		f[1][j][1]=ft[im][j-1][1];
		f[1][j][2]=ft[im][j  ][2];
		f[1][j][3]=ft[im][j+1][3];

		f[im][j][6]=ft[1][j-1][6];
		f[im][j][7]=ft[1][j  ][7];
		f[im][j][8]=ft[1][j+1][8];
	}
}
