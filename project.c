#include "poisson.h"
#include "functions.h"

int main(int argc, char *argv[]){

    int resolution = atoi(argv[1]);
    double Re = (double) atoi(argv[2]);

    // initialization of the program
    PetscInitialize(&argc, &argv, 0, 0);
    Sim_data *data = malloc(sizeof(Sim_data));
    init_sim_data(data, resolution, Re);
    //Poisson_data *data = malloc(sizeof(Poisson_data));
    //initialize_poisson_solver(data, resolution);

    int M = data->M;
    int N = data->N;
    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N; j++) {
            data->u[1+j + i*(N+2)] = j*0.1;
        }
    }

    set_boundary(data);
    write_fields(data, "data/test.txt");

    // End of the program
    //free_poisson_solver(data);
    //free(data);
    PetscFinalize();
    free_sim_data(data);
    free(data);
}
