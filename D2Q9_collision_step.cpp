#include "lbm_D2Q9.h"


void D2Q9_collision_step(double ***f, double ***ft, double **rho, double **u, double **v, double **uu, double **uv, double **vu, double **vv, double **type)
{
	#pragma omp parallel for
	for(int i=1; i<im+1; i++){
		for(int j=1; j<jm+1; j++){

			if(type[i][j]==FIELD)
			{
				double udotu=u[i][j]*u[i][j]+v[i][j]*v[i][j];

				double PI00=uu[i][j];
				double PI01=uv[i][j];
				double PI10=vu[i][j];
				double PI11=vv[i][j];

				double PIeq00=rho[i][j]*(u[i][j]*u[i][j]+1.0/3.0);
				double PIeq01=rho[i][j]*(u[i][j]*v[i][j]);
				double PIeq10=rho[i][j]*(v[i][j]*u[i][j]);
				double PIeq11=rho[i][j]*(v[i][j]*v[i][j]+1.0/3.0);

				double PIneq00=PI00-PIeq00;
				double PIneq01=PI01-PIeq01;
				double PIneq10=PI10-PIeq10;
				double PIneq11=PI11-PIeq11;

				for(int q=0; q<9; q++)
				{
					double edotu=ex[q]*u[i][j]+ey[q]*v[i][j];

					double feq=a[q]*rho[i][j]*(1.0+3.0*edotu+4.5*edotu*edotu-1.5*udotu);

					double Q00=(ex[q]*ex[q]-1.0/3.0);
					double Q01=(ex[q]*ey[q]);
					double Q10=(ey[q]*ex[q]);
					double Q11=(ey[q]*ey[q]-1.0/3.0);

					//double Quu=  Q00*u[i][j]*u[i][j]
					//			+Q01*u[i][j]*v[i][j]
					//			+Q10*v[i][j]*u[i][j]
					//			+Q11*v[i][j]*v[i][j];

					//double feq=a[q]*rho[i][j]*(1.0+3.0*edotu+4.5*Quu);

					double f1= 4.5*a[q]*(Q00*PIneq00
										+Q01*PIneq01
										+Q10*PIneq10
										+Q11*PIneq11);

					ft[i][j][q]=feq+(1.0-omega)*f1; // RBGK

					//ft[i][j][q]=(1.0-omega)*f[i][j][q]+omega*feq; // BGK
				}
			}
		}
	}
}