#include "poisson.h"
#include "functions.h"
#include "cases.h"

int main(int argc, char *argv[]) {

    int resolution = atoi(argv[1]);
    double Re = 500.0;
    double t_star_f = 50.0;
    int nb_frames = 400;
    char *folder = "data/video_test";

    PetscInitialize(&argc, &argv, 0, 0);

    case_oscillations(resolution, Re, t_star_f, nb_frames, folder);

    PetscFinalize();
}
