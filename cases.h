#ifndef _CASES_H_
#define _CASES_H_

void case_static(int resolution, double Re, double t_star_f, int nb_frames, char *folder);
void case_initial_perturbation(int resolution, double Re, double t_star_f, int nb_frames, char *folder);
void case_oscillations(int resolution, double Re, double t_star_f, int nb_frames, char *folder);
void case_oscillations_pert(int resolution, double Re, double t_star_f, int nb_frames, char *folder);

#endif
