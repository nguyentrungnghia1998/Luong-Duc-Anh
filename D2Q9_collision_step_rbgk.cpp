//------------------------------------------------------------------------------------------------
// Regularized Bhatnagar-Gross-Krook (RBGK) collision model
//------------------------------------------------------------------------------------------------


#include "lbm_D2Q9.h"


void D2Q9_collision_step_rbgk(double ***f, double ***fpc, double **type)
{
	#pragma omp parallel for
	for(int i=1; i<im+1; i++){
		for(int j=1; j<jm+1; j++){

			if(type[i][j]==FIELD)
			{
				double M00=0.0;
				double M10=0.0;
				double M01=0.0;
				double M11=0.0;
				double M20=0.0;
				double M02=0.0;

				for(int q=0; q<9; q++)
				{
					M00+=f[i][j][q];
					M10+=f[i][j][q]*ex[q];
					M01+=f[i][j][q]*ey[q];
					M11+=f[i][j][q]*ex[q]*ey[q];
					M20+=f[i][j][q]*ex[q]*ex[q];
					M02+=f[i][j][q]*ey[q]*ey[q];
				}

				double rho=M00;
				double ux=M10/rho;
				double uy=M01/rho;
				double PI00=M20;
				double PI01=M11;
				double PI10=M11;
				double PI11=M02;


				double feq[9];
				double uu=ux*ux+uy*uy;

				for(int q=0; q<9; q++)
				{
					double eu=ex[q]*ux+ey[q]*uy;
					feq[q]=a[q]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uu);
				}


				double f1[9];

				double PIeq00=rho*(ux*ux+1.0/3.0);
				double PIeq01=rho*(ux*uy);
				double PIeq10=rho*(uy*ux);
				double PIeq11=rho*(uy*uy+1.0/3.0);

				double PIneq00=PI00-PIeq00;
				double PIneq01=PI01-PIeq01;
				double PIneq10=PI10-PIeq10;
				double PIneq11=PI11-PIeq11;

				for(int q=0; q<9; q++)
				{
					double Q00=(ex[q]*ex[q]-1.0/3.0);
					double Q01=(ex[q]*ey[q]);
					double Q10=(ey[q]*ex[q]);
					double Q11=(ey[q]*ey[q]-1.0/3.0);

					f1[q]= 4.5*a[q]*(Q00*PIneq00+Q01*PIneq01
									+Q10*PIneq10+Q11*PIneq11);
				}


				for(int q=0; q<9; q++) fpc[i][j][q]=feq[q]+(1.0-omega)*f1[q];

			}
		}
	}
}