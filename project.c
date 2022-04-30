#include "poisson.h"
#include "functions.h"

int main(int argc, char *argv[]) {

    int resolution = atoi(argv[1]);
    double Re = (double) atoi(argv[2]);
    int nIt = 1;

    // initialization of the program
    PetscInitialize(&argc, &argv, 0, 0);
    Sim_data *sdata = malloc(sizeof(Sim_data));
    init_sim_data(sdata, resolution, Re);
    Poisson_data *pdata = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(pdata, resolution);

    int M = sdata->M;
    int N = sdata->N;

    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N; j++) {
            sdata->u[ j+1 + (N+2)*i ]     = 500*cos(2*M_PI*i/M);
            sdata->u_pre[ j+1 + (N+2)*i ] = 500*cos(2*M_PI*i/M);
        }
    }

    for (int i = 0; i < nIt; i++) {
        //compute_star(sdata);
        //set_boundary(sdata);
        //write_fields(sdata, "data/file.txt", 1);
        //switch_n(sdata);
        //mass_flow_condition(sdata);
        //poisson_solver(pdata, sdata);
        //update_pressure(sdata);
    }

    //write_fields(sdata, "data/file.txt", 0);

    // End of the program
    free_poisson_solver(pdata);
    free(pdata);
    PetscFinalize();
    free_sim_data(sdata);
    free(sdata);
}
