#include "lbm_D2Q9.h"


void D2Q9_macro_quantity(double ***f, double **rho, double **u, double **v, double **type)
{
	#pragma omp parallel for
	for(int i=1; i<im+1; i++){
		for(int j=1; j<jm+1; j++){

			if(type[i][j]==FIELD)
			{
				// raw moments
				double M00=0.0; // 0th order
				double M10=0.0; // 1st order
				double M01=0.0;
				//double M11=0.0; // 2nd order
				//double M20=0.0;
				//double M02=0.0;
				//double M21=0.0; // 3rd order
				//double M12=0.0;
				//double M22=0.0; // 4th order

				for(int q=0; q<9; q++)
				{
					M00+=f[i][j][q];
					M10+=f[i][j][q]*ex[q];
					M01+=f[i][j][q]*ey[q];
					//M11+=f[i][j][q]*ex[q]*ey[q];
					//M20+=f[i][j][q]*ex[q]*ex[q];
					//M02+=f[i][j][q]*ey[q]*ey[q];
					//M21+=f[i][j][q]*ex[q]*ex[q]*ey[q];
					//M12+=f[i][j][q]*ex[q]*ey[q]*ey[q];
					//M22+=f[i][j][q]*ex[q]*ex[q]*ey[q]*ey[q];
				}

				rho[i][j]=M00;
				u[i][j]=M10/rho[i][j];
				v[i][j]=M01/rho[i][j];
			}
			else /*if(type[i][j]==BOUND)*/
			{
				rho[i][j]=1.0;
				u[i][j]=0.0;
				v[i][j]=0.0;
			}
		}
	}
}