#include "poisson.h"
#include "functions.h"
#include "cases.h"

int main(int argc, char *argv[]) {

    int resolution = 50;
    double Re = 500.0;
    double t_star_f = 50.0;
    int nb_frames = 800;
    char *folder;

    int i = atoi(argv[1]);

    PetscInitialize(&argc, &argv, 0, 0);

    if (i == 1) {
        folder = "data/case1";
        case_static(resolution, Re, t_star_f, nb_frames, folder);
    }

    if (i == 2) {
        folder = "data/case2";
        case_initial_perturbation(resolution, Re, t_star_f, nb_frames, folder);
    }

    if (i == 3) {
        folder = "data/case3";
        case_oscillations(resolution, Re, t_star_f, nb_frames, folder);
    }

    if (i == 4) {
        folder = "data/case4";
        case_oscillations_pert(resolution, Re, t_star_f, nb_frames, folder);
    }

    PetscFinalize();
}
