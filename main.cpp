//----------------------------------------------------------------
// Lattice Boltzmann method (LBM) - D2Q9 lattice model
//----------------------------------------------------------------

//----------------------------------------------------------------
// Author: Kazuaki Uchibori
// E-mail: kazuakiuchibori@gmail.com
// (c) 2020 Kazuaki Uchibori
//----------------------------------------------------------------


#include "lbm_D2Q9.h"


int main(void)
{
	clock_t ttst,tted;    // elapsed time

	unsigned int n=0;    // time step
	unsigned int nout=0; // field output count
	double t=0.0;        // simulation time

	double ***f;         // particle distribution function (PDF)
	double ***fpc;       // PDF at post-collision state
	double **rho;        // macroscopic density
	double **u;          // macroscopic velocity x
	double **v;          // macroscopic velocity y
	double **type;       // celltype (FIELD or BOUND)

	D2Q9_alloc_pdf(f);
	D2Q9_alloc_pdf(fpc);
	D2Q9_alloc_qty(rho);
	D2Q9_alloc_qty(u);
	D2Q9_alloc_qty(v);
	D2Q9_alloc_qty(type);
	
	D2Q9_init_celltype(type);//sang
	D2Q9_init_field(rho,u,v,type);//thế
	D2Q9_init_pdf(f,fpc,rho,u,v);//sang
	D2Q9_macro_quantity(f,rho,u,v,type);//thế

	nout=D2Q9_output_vtk(nout,rho,u,v,type);

	printf("omega = %f (0 ~ 2) \n",omega);
	printf("Re = %f \n",re);

	while(t<tend)
	{
		ttst=clock();
		n++; printf("[%d] ",n); fflush(stdout);
		t+=dt; printf("t=%f ",t); fflush(stdout);

		D2Q9_collision_step_bgk(f,fpc,type);//Sơn
		//D2Q9_collision_step_rbgk(f,fpc,type);
		// D2Q9_collision_step_emrt(f,fpc,type);
		D2Q9_streaming_step(f,fpc,type);

		//D2Q9_domain_bc_periodic(f,fpc);
		D2Q9_domain_bc_slip(f,fpc);

		if(t>=nout*tout){
			D2Q9_macro_quantity(f,rho,u,v,type);
			nout=D2Q9_output_vtk(nout,rho,u,v,type);
			printf("u=%f ",u[im/2][jm/2]); fflush(stdout);
		}

		tted=clock();
		printf("%.3f sec\n",(double)(tted-ttst)/CLOCKS_PER_SEC); fflush(stdout);
		if((fopen("stop","r"))!=NULL) {
            break;
        }
	}
	
	D2Q9_free_pdf(f);
	D2Q9_free_pdf(fpc);
	D2Q9_free_qty(rho);
	D2Q9_free_qty(u);
	D2Q9_free_qty(v);
	D2Q9_free_qty(type);

	return 0;
}