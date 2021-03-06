#include <mpi.h>
#include <math.h>
#include "poisson.h"
#include "functions.h"

/*Modification to do :*/
/*    -Impose zero mass flow here by changing value of U_star*/
/*    -Fill vector rhs*/
void computeRHS(double *rhs, Sim_data *data) {

    mass_flow_condition(data);

    int N = data->N;
    int M = data->M;
    int idx;

    for (int i = 1; i < M+1; i++) {
        for (int j = 1; j < N+1; j++) {
            idx = (i-1)*N + j-1;
            rhs[idx] = divergence(data, i, j);
        }
    }

    rhs[0] = 0.0;
}

/*To call at each time step after computation of U_star. This function solves the poisson equation*/
/*and copies the solution of the equation into your vector Phi*/
/*Modification to do :*/
/*    - Change the call to computeRHS as you have to modify its prototype too*/
/*    - Copy solution of the equation into your vector PHI*/
void poisson_solver(Poisson_data *pdata, Sim_data *sdata) {

    /* Solve the linear system Ax = b for a 2-D poisson equation on a structured grid */
    int its;
    PetscInt rowStart, rowEnd;
    PetscScalar *rhs, *sol;

    KSP sles = pdata->sles;
    Vec b = pdata->b;
    Vec x = pdata->x;

    /* Fill the right-hand-side vector : b */
    VecGetOwnershipRange(b, &rowStart, &rowEnd);
    VecGetArray(b, &rhs);
    computeRHS(rhs, sdata); /*MODIFY THE PROTOTYPE HERE*/
    VecRestoreArray(b, &rhs);

    /*Solve the linear system of equations */
    KSPSolve(sles, b, x);
    KSPGetIterationNumber(sles, &its);
    //PetscPrintf(PETSC_COMM_WORLD, "Solution to Poisson eqn in %d iterations \n", its);

    VecGetArray(x, &sol);

    int r;
    for (r = rowStart; r < rowEnd; r++){
        sdata->phi[r] = sol[r];
    }

    VecRestoreArray(x, &sol);
    free(sol);
    free(rhs);
}

/*This function is called only once during the simulation, i.e. in initialize_poisson_solver.*/
/*In its current state, it inserts unity on the main diagonal.*/
/*More than probably, you should need to add arguments to the prototype ... .*/
/*Modification to do in this function : */
/*   -Insert the correct factor in matrix A*/
void computeLaplacianMatrix(Mat A, int resolution, double h, double dt) {

    int N = 5*resolution;
    int M = 15*resolution;
    int idx;

    // Interior points
    for (int i = 1; i < M-1; i++) {
        for (int j = 1; j < N-1; j++) {
            idx = j + i*N;
            MatSetValue(A, idx, idx, -4.0*dt/h, INSERT_VALUES);
            MatSetValue(A, idx, idx+1, 1.0*dt/h, INSERT_VALUES);
            MatSetValue(A, idx, idx-1, 1.0*dt/h, INSERT_VALUES);
            MatSetValue(A, idx, idx+N, 1.0*dt/h, INSERT_VALUES);
            MatSetValue(A, idx, idx-N, 1.0*dt/h, INSERT_VALUES);
        }
    }

    // Boundary points
    for (int j = 1; j < N-1; j++) {
        // left boundary
        idx = j;
        MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx-1, 1.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx+1, 1.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx+N, 1.0*dt/h, INSERT_VALUES);

        // right boundary
        idx = j + (M-1)*N;
        MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx-1, 1.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx+1, 1.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx-N, 1.0*dt/h, INSERT_VALUES);
    }

    for (int i = 1; i < M-1; i++) {
        // lower boundary
        idx = i*N;
        MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx+1, 1.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx-N, 1.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx+N, 1.0*dt/h, INSERT_VALUES);

        // upper boundary
        idx = i*N + N-1;
        MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx-1, 1.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx-N, 1.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx+N, 1.0*dt/h, INSERT_VALUES);
    }

    // Corners
    MatSetValue(A, 0, 0, 1.0, INSERT_VALUES);
    idx = N-1;
    MatSetValue(A, idx, idx, -2.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx-1, 1.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx+N, 1.0*dt/h, INSERT_VALUES);
    idx = (M-1)*N;
    MatSetValue(A, idx, idx, -2.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx+1, 1.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx-N, 1.0*dt/h, INSERT_VALUES);
    idx = N*M-1;
    MatSetValue(A, idx, idx, -2.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx-1, 1.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx-N, 1.0*dt/h, INSERT_VALUES);

    // Autour de l'objet
    int i_start = 3*resolution;
    int i_to    = 8*resolution;
    int j_start = 2*resolution;
    int j_to    = 3*resolution;

    for (int i = i_start; i < i_to; i++) {
        // lower side
        idx = j_start-1 + i*N;
        MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx+1, 0.0, INSERT_VALUES);
        // upper side
        idx = j_to + i*N;
        MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx-1, 0.0, INSERT_VALUES);
    }
    for (int j = j_start; j < j_to; j++) {
        // left side
        idx = j + (i_start-1)*N;
        MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx+N, 0.0, INSERT_VALUES);
        // right side
        idx = j + i_to*N;
        MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
        MatSetValue(A, idx, idx-N, 0.0, INSERT_VALUES);
    }
    // coins du bloc
    idx = (resolution*3 - 1)*N + (2*resolution - 1);    // bas gauche
    MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx+1, 0.0, INSERT_VALUES);
    MatSetValue(A, idx, idx+N, 0.0, INSERT_VALUES);
    idx = (resolution*3 - 1)*N + (3*resolution);        // haut gauche
    MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx-1, 0.0, INSERT_VALUES);
    MatSetValue(A, idx, idx+N, 0.0, INSERT_VALUES);
    idx = (resolution*8)*N + (2*resolution - 1);        // bas droite
    MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx+1, 0.0, INSERT_VALUES);
    MatSetValue(A, idx, idx-N, 0.0, INSERT_VALUES);
    idx = (resolution*8)*N + (3*resolution);            // haut droite
    MatSetValue(A, idx, idx, -3.0*dt/h, INSERT_VALUES);
    MatSetValue(A, idx, idx-1, 0.0, INSERT_VALUES);
    MatSetValue(A, idx, idx-N, 0.0, INSERT_VALUES);
}

/*To call during the initialization of your solver, before the begin of the time loop*/
/*Maybe you should need to add an argument to specify the number of unknows*/
/*Modification to do in this function :*/
/*   -Specify the number of unknows*/
/*   -Specify the number of non-zero diagonals in the sparse matrix*/
PetscErrorCode initialize_poisson_solver(Poisson_data* data, int resolution, double h, double dt) {
    PetscInt rowStart; /*rowStart = 0*/
    PetscInt rowEnd; /*rowEnd = the number of unknows*/
    PetscErrorCode ierr;

    int N = 5*resolution;
    int M = 15*resolution;
    int nphi = N*M;

    /* Create the right-hand-side vector : b */
    VecCreate(PETSC_COMM_WORLD, &(data->b));
    VecSetSizes(data->b, PETSC_DECIDE, nphi);
    VecSetType(data->b,VECSTANDARD);

    /* Create the solution vector : x */
    VecCreate(PETSC_COMM_WORLD, &(data->x));
    VecSetSizes(data->x, PETSC_DECIDE, nphi);
    VecSetType(data->x,VECSTANDARD);

    /* Create and assemble the Laplacian matrix : A  */
    MatCreate(PETSC_COMM_WORLD, &(data->A));
    MatSetSizes(data->A, PETSC_DECIDE, PETSC_DECIDE, nphi , nphi);
    MatSetType(data->A, MATAIJ);
    MatSeqAIJSetPreallocation(data->A, 5, NULL);
    MatGetOwnershipRange(data->A, &rowStart, &rowEnd);

    computeLaplacianMatrix(data->A, resolution, h, dt);
    ierr = MatAssemblyBegin(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(data->A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /* Create the Krylov context */
    KSPCreate(PETSC_COMM_WORLD, &(data->sles));
    KSPSetOperators(data->sles, data->A, data->A);
    KSPSetType(data->sles,KSPGMRES); //KSPGMRES seems the best, it will not be used if PC LU.
    PC prec;
    KSPGetPC(data->sles, &prec);
    PCSetType(prec,PCLU);
    KSPSetFromOptions(data->sles); // to uncomment if we want to specify the solver to use in command line. Ex: mpirun -ksp_type gmres -pc_type gamg
    KSPSetTolerances(data->sles, 1.e-12, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetReusePreconditioner(data->sles,PETSC_TRUE);
    KSPSetUseFischerGuess(data->sles,1,4);
    KSPGMRESSetPreAllocateVectors(data->sles);

    PetscPrintf(PETSC_COMM_WORLD, "Assembly of Matrix and Vectors is done \n");

    return ierr;
}

/*To call after the simulation to free the vectors needed for Poisson equation*/
/*Modification to do : nothing */
void free_poisson_solver(Poisson_data* data){
    MatDestroy(&(data->A));
    VecDestroy(&(data->b));
    VecDestroy(&(data->x));
    KSPDestroy(&(data->sles));
}
