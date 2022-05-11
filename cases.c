#include <math.h>
#include <time.h>
#include "poisson.h"
#include "functions.h"
#include "cases.h"

void case_static(int resolution, double Re, double t_star_f, int nb_frames, char *folder) {

    // initialization of the program
    Sim_data *sdata = malloc(sizeof(Sim_data));
    init_sim_data(sdata, resolution, Re);
    Poisson_data *pdata = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(pdata, resolution, sdata->h, sdata->dt);

    int nIt = (int) (t_star_f * sdata->H_box / sdata->U_inf / sdata->dt);

    printf("# of iterations = %d\n", nIt);
    printf("dt = %fs\n", sdata->dt);

    int M = sdata->M;
    int N = sdata->N;

    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N+2; j++) {
            sdata->u[ j + (N+2)*i ]     = sdata->U_inf;
            sdata->u_pre[ j + (N+2)*i ] = sdata->U_inf;
        }
    }
    set_boundary(sdata);

    int step = nIt / nb_frames;
    char *filename = malloc(50*sizeof(char));
    clock_t t = clock();

    for (int i = 0; i < nIt; i++) {

        sdata->t_star += sdata->dt * sdata->U_inf / sdata->H_box;

        compute_star(sdata);            // calcul u* à l'intérieur
        switch_n(sdata);                // met u_n dans u_n-1
        set_boundary(sdata);            // calcul u* sur la frontière et les points fantomes
        poisson_solver(pdata, sdata);   // calcul phi
        compute_next(sdata);            // calcul u_n+1 (et met u* sur la frontière)
        update_pressure(sdata);         // met à jour la pression

        if (i%100 == 0) {
            printf("iteration %d\n", i);
            printf("t* = %.2f\n", sdata->t_star);
            printf("time remaining: %.1f sec\n", time_remaining(i, nIt, &t));
        }
        if (i%step == 0) {
            sprintf(filename, "%s/step_%d.txt", folder, i/step + 1);
            write_fields(sdata, filename);
        }
    }

    // End of the program
    free_poisson_solver(pdata);
    free(pdata);
    free_sim_data(sdata);
    free(sdata);
}

void case_initial_perturbation(int resolution, double Re, double t_star_f, int nb_frames, char *folder) {

    // initialization of the program
    Sim_data *sdata = malloc(sizeof(Sim_data));
    init_sim_data(sdata, resolution, Re);
    Poisson_data *pdata = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(pdata, resolution, sdata->h, sdata->dt);

    int nIt = (int) (t_star_f * sdata->H_box / sdata->U_inf / sdata->dt);

    printf("# of iterations = %d\n", nIt);
    printf("dt = %fs\n", sdata->dt);

    int M = sdata->M;
    int N = sdata->N;

    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N+2; j++) {
            sdata->u[ j + (N+2)*i ]     = sdata->U_inf;
            sdata->u_pre[ j + (N+2)*i ] = sdata->U_inf;
        }
    }
    set_boundary(sdata);

    int step = nIt / nb_frames;
    char *filename = malloc(50*sizeof(char));
    clock_t t = clock();

    for (int i = 0; i < nIt; i++) {

        sdata->t_star += sdata->dt * sdata->U_inf / sdata->H_box;

        if (sdata->t_star <= 1.0) {
            sdata->v_mesh = sdata->U_inf * sin(2.0 * M_PI * sdata->t_star) / 10.0;
            sdata->y_mesh = (1.0 - cos(2.0 * M_PI * sdata->t_star)) / (20.0 * M_PI);
        } else {
            sdata->v_mesh = 0.0;
            sdata->y_mesh = 0.0;
        }

        compute_star(sdata);
        switch_n(sdata);
        set_boundary(sdata);
        poisson_solver(pdata, sdata);   // calcul phi
        compute_next(sdata);            // calcul u_n+1 (et met u* sur la frontière)
        update_pressure(sdata);         // met à jour la pression

        if (i%100 == 0) {
            printf("iteration %d\n", i);
            printf("t* = %.2f\n", sdata->t_star);
            printf("time remaining: %.1f sec\n", time_remaining(i, nIt, &t));
        }
        if (i%step == 0) {
            sprintf(filename, "%s/step_%d.txt", folder, i/step + 1);
            write_fields(sdata, filename);
        }

        sdata->v_mesh_prev = sdata->v_mesh;
    }

    // End of the program
    free_poisson_solver(pdata);
    free(pdata);
    free_sim_data(sdata);
    free(sdata);
}

void case_oscillations(int resolution, double Re, double t_star_f, int nb_frames, char *folder) {

    // initialization of the program
    Sim_data *sdata = malloc(sizeof(Sim_data));
    init_sim_data(sdata, resolution, Re);
    Poisson_data *pdata = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(pdata, resolution, sdata->h, sdata->dt);

    int nIt = (int) (t_star_f * sdata->H_box / sdata->U_inf / sdata->dt);

    printf("# of iterations = %d\n", nIt);
    printf("dt = %fs\n", sdata->dt);

    int M = sdata->M;
    int N = sdata->N;

    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N+2; j++) {
            sdata->u[ j + (N+2)*i ]     = sdata->U_inf;
            sdata->u_pre[ j + (N+2)*i ] = sdata->U_inf;
        }
    }
    set_boundary(sdata);

    int step = nIt / nb_frames;
    char *filename = malloc(50*sizeof(char));
    clock_t t = clock();

    for (int i = 0; i < nIt; i++) {

        sdata->t_star += sdata->dt * sdata->U_inf / sdata->H_box;
        sdata->u_mesh = sdata->U_inf * sdata->alpha * sin(2.0 * M_PI * sdata->St * sdata->t_star);
        sdata->x_mesh = sdata->alpha / (2.0 * M_PI * sdata->St) * (1.0 - cos(2.0 * M_PI * sdata->St * sdata->t_star));

        compute_star(sdata);
        switch_n(sdata);                // met u_n dans u_n-1
        set_boundary(sdata);            // calcul u* sur la frontière et les points fantomes
        poisson_solver(pdata, sdata);   // calcul phi
        compute_next(sdata);            // calcul u_n+1 (et met u* sur la frontière)
        update_pressure(sdata);         // met à jour la pression

        if (i%100 == 0) {
            printf("iteration %d\n", i);
            printf("t* = %.2f\n", sdata->t_star);
            printf("time remaining: %.1f sec\n", time_remaining(i, nIt, &t));
        }
        if (i%step == 0) {
            sprintf(filename, "%s/step_%d.txt", folder, i/step + 1);
            write_fields(sdata, filename);
        }

        sdata->u_mesh_prev = sdata->u_mesh;
    }

    // End of the program
    free_poisson_solver(pdata);
    free(pdata);
    free_sim_data(sdata);
    free(sdata);
}

void case_oscillations_pert(int resolution, double Re, double t_star_f, int nb_frames, char *folder) {

    // initialization of the program
    Sim_data *sdata = malloc(sizeof(Sim_data));
    init_sim_data(sdata, resolution, Re);
    Poisson_data *pdata = malloc(sizeof(Poisson_data));
    initialize_poisson_solver(pdata, resolution, sdata->h, sdata->dt);

    int nIt = (int) (t_star_f * sdata->H_box / sdata->U_inf / sdata->dt);

    printf("# of iterations = %d\n", nIt);
    printf("dt = %fs\n", sdata->dt);

    int M = sdata->M;
    int N = sdata->N;

    for (int i = 0; i < M+1; i++) {
        for (int j = 0; j < N+2; j++) {
            sdata->u[ j + (N+2)*i ]     = sdata->U_inf;
            sdata->u_pre[ j + (N+2)*i ] = sdata->U_inf;
        }
    }
    set_boundary(sdata);

    int step = nIt / nb_frames;
    char *filename = malloc(50*sizeof(char));
    clock_t t = clock();

    for (int i = 0; i < nIt; i++) {

        sdata->t_star += sdata->dt * sdata->U_inf / sdata->H_box;
        sdata->u_mesh = sdata->U_inf * sdata->alpha * sin(2.0 * M_PI * sdata->St * sdata->t_star);
        sdata->x_mesh = sdata->alpha / (2.0 * M_PI * sdata->St) * (1.0 - cos(2.0 * M_PI * sdata->St * sdata->t_star));

        if (sdata->t_star <= 1.0) {
            sdata->v_mesh = sdata->U_inf * sin(2.0 * M_PI * sdata->t_star) / 10.0;
            sdata->y_mesh = (1.0 - cos(2.0 * M_PI * sdata->t_star)) / (20.0 * M_PI);
        } else {
            sdata->v_mesh = 0.0;
            sdata->y_mesh = 0.0;
        }

        compute_star(sdata);
        switch_n(sdata);
        set_boundary(sdata);
        poisson_solver(pdata, sdata);   // calcul phi
        compute_next(sdata);            // calcul u_n+1 (et met u* sur la frontière)
        update_pressure(sdata);         // met à jour la pression

        if (i%100 == 0) {
            printf("iteration %d\n", i);
            printf("t* = %.2f\n", sdata->t_star);
            printf("time remaining: %.1f sec\n", time_remaining(i, nIt, &t));
        }
        if (i%step == 0) {
            sprintf(filename, "%s/step_%d.txt", folder, i/step + 1);
            write_fields(sdata, filename);
        }

        sdata->u_mesh_prev = sdata->u_mesh;
        sdata->v_mesh_prev = sdata->v_mesh;
    }

    // End of the program
    free_poisson_solver(pdata);
    free(pdata);
    free_sim_data(sdata);
    free(sdata);
}
