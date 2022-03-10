#include "lbm_D2Q9.h"


void D2Q9_alloc_pdf(double ***&f)
{
	f=(double***)malloc(sizeof(double**)*(im+2));
	for(int i=0;i<(im+2);i++){
		f[i]=(double**)malloc(sizeof(double*)*(jm+2));
		for(int j=0;j<(jm+2);j++){
			f[i][j]=(double*)malloc(sizeof(double)*9);
		}
	}
}

void D2Q9_alloc_qty(double **&f)
{
	f=(double**)malloc(sizeof(double*)*(im+2));
	for(int i=0;i<(im+2);i++){
		f[i]=(double*)malloc(sizeof(double)*(jm+2));
	}
}

void D2Q9_free_pdf(double ***f)
{
	for(int i=0;i<(im+2);i++){
		for(int j=0;j<(jm+2);j++){
			free(f[i][j]);
		}
		free(f[i]);
	}
	free(f);
}

void D2Q9_free_qty(double **f)
{
	for(int i=0;i<(im+2);i++){
		free(f[i]);
	}
	free(f);
}