//------------------------------------------------------------------------------------------------
// Entropic multi-relaxation time (EMRT) collision model
// Karlin-Bosch-Chikatamarla (KBC) collision model
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
// Reference
//
// 1. F. Bosch, S. S. Chikatamarla and I. Karlin, "Entropic multi-relaxation lattice Boltzmann
//    scheme for turbulent flows", Physical Review E 92 (4), 043309, (2015).
//
// 2. C. Coreixas, B. Chopard and J. Latt, "Comprehensive comparison of collision models in the 
//    lattice Boltzmann framework: Theoretical investigations", Physical Review E 100, 033305, 
//    (2019).
//------------------------------------------------------------------------------------------------


#include "lbm_D2Q9.h"


void D2Q9_collision_step_emrt(double ***f, double ***fpc, double **type)
{
	double beta=omega*0.5;
	double gamma0=1.0/beta;

	#pragma omp parallel for
	for(int i=1; i<im+1; i++){
		for(int j=1; j<jm+1; j++){

			if(type[i][j]==FIELD)
			{
				// raw moment
				double rhoM00=0.0;
				double rhoM10=0.0;
				double rhoM01=0.0;
				double rhoM11=0.0;
				double rhoM20=0.0;
				double rhoM02=0.0;

				for(int q=0; q<9; q++)
				{
					rhoM00+=f[i][j][q];
					rhoM10+=f[i][j][q]*ex[q];
					rhoM01+=f[i][j][q]*ey[q];
					rhoM11+=f[i][j][q]*ex[q]*ey[q];
					rhoM20+=f[i][j][q]*ex[q]*ex[q];
					rhoM02+=f[i][j][q]*ey[q]*ey[q];
				}

				double rho=rhoM00;


				// shear part: s
				double s[9];
				s[0]=-rhoM20-rhoM02;
				s[1]= rhoM11*0.25;
				s[2]= rhoM20*0.50;
				s[3]=-rhoM11*0.25;
				s[4]= rhoM02*0.50;
				s[5]= rhoM02*0.50;
				s[6]=-rhoM11*0.25;
				s[7]= rhoM20*0.50;
				s[8]= rhoM11*0.25;


				// equilibrium: feq
				double feq[9];
				double ux=rhoM10/rho;
				double uy=rhoM01/rho;
				double uu=ux*ux+uy*uy;

				for(int q=0; q<9; q++)
				{
					double eu=ex[q]*ux+ey[q]*uy;
					feq[q]=a[q]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uu);
				}


				// raw moment w/ equilibrium
				double rhoMeq00=0.0;
				double rhoMeq11=0.0;
				double rhoMeq20=0.0;
				double rhoMeq02=0.0;

				for(int q=0; q<9; q++)
				{
					rhoMeq00+=feq[q];
					rhoMeq11+=feq[q]*ex[q]*ey[q];
					rhoMeq20+=feq[q]*ex[q]*ex[q];
					rhoMeq02+=feq[q]*ey[q]*ey[q];
				}


				// equilibrium of shear part: seq
				double seq[9];
				seq[0]=-rhoMeq20-rhoMeq02;
				seq[1]= rhoMeq11*0.25;
				seq[2]= rhoMeq20*0.50;
				seq[3]=-rhoMeq11*0.25;
				seq[4]= rhoMeq02*0.50;
				seq[5]= rhoMeq02*0.50;
				seq[6]=-rhoMeq11*0.25;
				seq[7]= rhoMeq20*0.50;
				seq[8]= rhoMeq11*0.25;


				// derivation of shear part: ds
				double ds[9];

				for(int q=0; q<9; q++) ds[q]=s[q]-seq[q];

				
				// derivation of higher-order part: dh
				double dh[9];

				for(int q=0; q<9; q++) dh[q]=f[i][j][q]-feq[q]-ds[q];


				// relaxation parameter: gamma
				double dsdh=0.0;
				double dhdh=0.0;

				for(int q=0; q<9; q++)
				{
					dsdh+=ds[q]*dh[q]/feq[q];
					dhdh+=dh[q]*dh[q]/feq[q];
				}

				double gamma=gamma0-(2.0-gamma0)*(dsdh/(dhdh+SMALL_VAL));


				// relaxation
				for(int q=0; q<9; q++) fpc[i][j][q]=f[i][j][q]-beta*(2.0*ds[q]+gamma*dh[q]);
			}
		}
	}
}