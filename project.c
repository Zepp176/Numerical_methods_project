#include "poisson.h"
#include "functions.h"

int main(int argc, char *argv[]) {

    int resolution = atoi(argv[1]);
    double Re = (double) atoi(argv[2]);

    // initialization of the program
    PetscInitialize(&argc, &argv, 0, 0);
    Sim_data *sdata = malloc(sizeof(Sim_data));
    init_sim_data(sdata, resolution, Re);
    Poisson_data *pdata = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(pdata, resolution, sdata->h, sdata->dt);

    int nIt = (int) atoi(argv[3]) * sdata->H_box / sdata->U_inf / sdata->dt;

    printf("dt = %f\n", sdata->dt);
    printf("h = %f\n", sdata->h);
    printf("t* = %f\n", nIt*sdata->dt*sdata->U_inf/sdata->H_box);
    printf("n = %d\n", nIt);

    int M = sdata->M;
    int N = sdata->N;

    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N+2; j++) {
            sdata->u[ j + (N+2)*i ]     = sdata->U_inf;
            sdata->u_pre[ j + (N+2)*i ] = sdata->U_inf;
        }
    }
    set_boundary(sdata);

    int nb_frames = 480;
    int step = nIt / nb_frames;
    char *str = malloc(50*sizeof(char));

    for (int i = 0; i < nIt; i++) {
        compute_star(sdata);            // calcul u* à l'intérieur
        switch_n(sdata);                // met u_n dans u_n-1
        set_boundary(sdata);            // calcul u* sur la frontière et les points fantomes
        poisson_solver(pdata, sdata);   // calcul phi
        compute_next(sdata);            // calcul u_n+1 (et met u* sur la frontière)
        update_pressure(sdata);         // met à jour la pression

        if (i%100 == 0) {
            printf("iteration %d\n", i);
        }
        if (i%step == 0) {
            sprintf(str, "data/video/step_%d.txt", i/step + 1);
            write_fields(sdata, str, 0);
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
