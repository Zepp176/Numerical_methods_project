#ifndef _POISSON_H_
#define _POISSON_H_

/*To include in the file in which you will call initialize_poisson_solver and poisson_solver*/

#include <petsc.h>
#include <petscsys.h>
#include "functions.h"

//Structure storing petsc vectors

typedef struct {

	Vec b;
	Vec x;
	Mat A;
	KSP sles;

} Poisson_data;

PetscErrorCode initialize_poisson_solver(Poisson_data* data, int resolution, double h, double dt);
void poisson_solver(Poisson_data *pdata, Sim_data *sdata);
void free_poisson_solver(Poisson_data* data);
void computeRHS(double *rhs, Sim_data *data);

#endif
