#ifndef LBM_D2Q9_H
#define LBM_D2Q9_H

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define SMALL_VAL 1.0E-20
#define LARGE_VAL 1.0E+20

#define FIELD 0.0
#define BOUND 1.0

static const int im=1000;
static const int jm=500;

static const double dr=1.0;
static const double dt=1.0;

static const double tend=dt*1e6;
static const double tout=dt*2e3;

static const double u0=5.0e-2;
static const double d0=dr*200;
static const double nu=2.0e-3;
static const double re=u0*d0/nu;

static const double cs=1.0/sqrt(3.0);

static const double nu_dash=nu*dt/(dr*dr);
static const double tau=3.0*nu_dash+0.5;
static const double omega=1.0/tau;

static const double ex[9]={ 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,-1.0,-1.0,-1.0 }; // lattice velocity x
static const double ey[9]={ 0.0, 1.0, 0.0,-1.0, 1.0,-1.0, 1.0, 0.0,-1.0 }; // lattice velocity y
static const double a[9]={4.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0 };

void D2Q9_alloc_pdf(double ***&f);
void D2Q9_alloc_qty(double **&f);
void D2Q9_init_celltype(double **type);
void D2Q9_init_field(double **rho, double **u, double **v, double **type);
void D2Q9_init_pdf(double ***f, double ***fpc, double **rho, double **u, double **v);
void D2Q9_collision_step_bgk(double ***f, double ***fpc, double **type);
void D2Q9_collision_step_rbgk(double ***f, double ***fpc, double **type);
void D2Q9_collision_step_emrt(double ***f, double ***fpc, double **type);
void D2Q9_streaming_step(double ***f, double ***fpc, double **type);
void D2Q9_macro_quantity(double ***f, double **rho, double **u, double **v, double **type);
void D2Q9_domain_bc_slip(double ***f, double ***ft);
void D2Q9_domain_bc_periodic(double ***f, double ***ft);
int D2Q9_output_vtk(int nout, double **rho, double **u, double **v, double **type);
void D2Q9_free_pdf(double ***f);
void D2Q9_free_qty(double **f);


#endif
