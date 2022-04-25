#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

typedef struct {

    double *u;
    double *v;
    double *u_pre;
    double *v_pre;
    double *u_star;
    double *v_star;
    double *P;
    double *phi;

    int M;
    int N;

    double h;
    double U_inf;
    double nu;
    double H_box;
    double dt;

} Sim_data;

double get_u(Sim_data *data, int t, float idx_i, float idx_j);
double get_v(Sim_data *data, int t, float idx_i, float idx_j);
double get_a(Sim_data *data, int t, float i, float j);
double get_b(Sim_data *data, int t, float i, float j);
double get_grad_P(Sim_data *data, float idx_i, float idx_j);
double get_laplacian(Sim_data *data, float idx_i, float idx_j);
void init_sim_data(Sim_data *data, int res, double Re);
void free_sim_data(Sim_data *data);
void write_fields(Sim_data *data, char *filename);
void set_boundary(Sim_data *data);
void compute_star(Sim_data *data);
void switch_n(Sim_data *data);
double divergence(Sim_data *data, int i, int j);
void update_pressure(Sim_data *data);

#endif
