//------------------------------------------------------------------------------------------------
// Bhatnagar-Gross-Krook (BGK) collision model
//------------------------------------------------------------------------------------------------


#include "lbm_D2Q9.h"


void D2Q9_collision_step_bgk(double ***f, double ***fpc, double **type)
{
	#pragma omp parallel for
	for(int i=1; i<im+1; i++){
		for(int j=1; j<jm+1; j++){

			if(type[i][j]==FIELD)
			{
				double M00=0.0;
				double M10=0.0;
				double M01=0.0;

				for(int q=0; q<9; q++)
				{
					M00+=f[i][j][q];
					M10+=f[i][j][q]*ex[q];
					M01+=f[i][j][q]*ey[q];
				}

				double feq[9];
				double rho=M00;
				double ux=M10/rho;
				double uy=M01/rho;
				double uu=ux*ux+uy*uy;

				for(int q=0; q<9; q++)
				{
					double eu=ex[q]*ux+ey[q]*uy;
					feq[q]=a[q]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uu);
				}

				for(int q=0; q<9; q++) fpc[i][j][q]=(1.0-omega)*f[i][j][q]+omega*feq[q];

			}
		}
	}
}