#include "lbm_D2Q9.h"


void D2Q9_init_field(double **rho, double **u, double **v, double **type)
{
	#pragma omp parallel for
	for(int i=0; i<im+2; i++){
		for(int j=0; j<jm+2; j++){

			rho[i][j]=1.0;
			u[i][j]=u0;
			v[i][j]=0.0;
		}
	}

	//#pragma omp parallel for
	//for(int i=0; i<im+2; i++){
	//	for(int j=0; j<jm+2; j++){

	//		if(j>jm/2) u[i][j]=u0;
	//		else u[i][j]=-u0;
	//		u[im/2][jm/2]=0.0;
	//	}
	//}

	#pragma omp parallel for
	for(int i=0; i<im+2; i++){
		for(int j=0; j<jm+2; j++){

			if(type[i][j]==BOUND)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
			}
		}
	}
}