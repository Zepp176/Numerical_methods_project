#include "poisson.h"
#include "functions.h"

int main(int argc, char *argv[]) {

    int resolution = atoi(argv[1]);
    double Re = (double) atoi(argv[2]);
    int nIt = 4000;

    // initialization of the program
    PetscInitialize(&argc, &argv, 0, 0);
    Sim_data *sdata = malloc(sizeof(Sim_data));
    init_sim_data(sdata, resolution, Re);
    Poisson_data *pdata = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(pdata, resolution, sdata->h, sdata->dt);

    printf("dt = %f\n", sdata->dt);
    printf("h = %f\n", sdata->h);

    int M = sdata->M;
    int N = sdata->N;

    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N; j++) {
            sdata->u[ j+1 + (N+2)*i ]     = sdata->U_inf*cos(2*M_PI*i/M);
            sdata->u_pre[ j+1 + (N+2)*i ] = sdata->U_inf*cos(2*M_PI*i/M);
        }
    }

    for (int i = 0; i < nIt; i++) {
        compute_star(sdata);
        set_boundary(sdata);
        mass_flow_condition(sdata);
        poisson_solver(pdata, sdata);
        switch_n(sdata);
        compute_next(sdata);
        update_pressure(sdata);

        if (i%100 == 0) {
            printf("iteration %d\n", i);
        }
    }

    write_fields(sdata, "data/file.txt", 0);

    // End of the program
    free_poisson_solver(pdata);
    free(pdata);
    PetscFinalize();
    free_sim_data(sdata);
    free(sdata);
}
