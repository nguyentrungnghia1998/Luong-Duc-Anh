#include "lbm_D2Q9.h"


void D2Q9_init_celltype(double **type)
{
	#pragma omp parallel for
	for(int i=0; i<im+2; i++){
		for(int j=0; j<jm+2; j++){
			type[i][j]=BOUND;
		}
	}

	#pragma omp parallel for
	for(int i=0; i<im+2; i++){
		for(int j=1; j<jm+1; j++){
			type[i][j]=FIELD;
		}
	}


	// flatplate
	//#pragma omp parallel for
	//for(int i=401; i<=401; i++){
	//	for(int j=421; j<=620; j++){
	//		type[i][j]=BOUND;
	//	}
	//}

	// circle
	double rad=dr*(jm/10);
	double cx=dr*(im/5);
	double cy=dr*(jm/2)+(jm/50);
	#pragma omp parallel for
	for(int i=0; i<im+2; i++){
		for(int j=0; j<jm+2; j++){

			double rx=dr*i-dr*0.5;
			double ry=dr*j-dr*0.5;

			if(sqrt((cx-rx)*(cx-rx)+(cy-ry)*(cy-ry))<=rad)
			{
				type[i][j]=BOUND;
			}
		}
	}
}
