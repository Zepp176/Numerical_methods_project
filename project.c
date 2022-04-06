#include "poisson.h"
#include "functions.h"

int main(int argc, char *argv[]){

    int resolution = atoi(argv[1]);
    int N = 5*resolution; int M = 15*resolution;

    double *u = malloc((N+2)*(M+1)*sizeof(double));
    double *v = malloc((N+1)*(M+2)*sizeof(double));
    double *v_star = malloc((N+1)*(M+2)*sizeof(double));
    double *u_pre = malloc((N+2)*(M+1)*sizeof(double));
    double *v_pre = malloc((N+1)*(M+2)*sizeof(double));
    double *P = malloc(N*M*sizeof(double));

    for (int i = 0; i < M+2; i++) {
        for (int j = 0; j < N+1; j++) {
            v[j+i*(N+1)] = j+i*(N+1);
        }
    }

    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N+2; j++) {
            u[j+i*(N+2)] = j+i*(N+2);
        }
    }

    free(u); free(v); free(v_star); free(u_pre); free(v_pre); free(P);

    /*
    // initialization of the program
    PetscInitialize(&argc, &argv, 0, 0);
    Poisson_data *data = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(data, resolution);
    PetscScalar *sol;
    FILE *f = fopen("data/file.txt", "w");
    if (f == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    poisson_solver(data, resolution);

    VecGetArray(data->x, &sol);

    if (f == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(f, "%d\n", resolution);

    for (int i = 0; i < resolution*5*resolution*15; i++) {
        fprintf(f, "%f\n", sol[i]);
    }

    // End of the program
    fclose(f);
    VecRestoreArray(data->x, &sol);
    free_poisson_solver(data);
    free(data);

    PetscFinalize();*/
}
