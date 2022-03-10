#include "lbm_D2Q9.h"


void D2Q9_init_pdf(double ***f, double ***fpc, double **rho, double **u, double **v)
{
	#pragma omp parallel for
	for(int i=0; i<im+2; i++){
		for(int j=0; j<jm+2; j++){

			double uu=u[i][j]*u[i][j]+v[i][j]*v[i][j];

			for(int q=0; q<9; q++)
			{
				double eu=ex[q]*u[i][j]+ey[q]*v[i][j];

				double feq=a[q]*(1.0+3.0*eu+4.5*eu*eu-1.5*uu);

				f[i][j][q]=feq;
				fpc[i][j][q]=feq;
			}
		}
	}

	#pragma omp parallel for
	for(int i=0; i<im+2; i++){
		for(int j=0; j<jm+2; j++){

			for(int q=0; q<9; q++)
			{
				fpc[i][j][q]=f[i][j][q];
			}
		}
	}
}