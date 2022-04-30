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
    int res;

} Sim_data;

double get_u(Sim_data *data, int t, float idx_i, float idx_j);
double get_v(Sim_data *data, int t, float idx_i, float idx_j);
double get_a(double *u, double *v, int M, int N, double h, int i, int j);
double get_b(double *u, double *v, int M, int N, double h, int i, int j);
double get_gradu(double *P, int M, int N, double h, int i, int j);
double get_gradv(double *P, int M, int N, double h, int i, int j);
double get_lapu(double *u, int M, int N, double h, int i, int j);
double get_lapv(double *v, int M, int N, double h, int i, int j);
double get_RHSu(Sim_data *data, int i, int j);
double get_RHSv(Sim_data *data, int i, int j);
void compute_star(Sim_data *data);
void compute_next(Sim_data *data);
void init_sim_data(Sim_data *data, int res, double Re);
void free_sim_data(Sim_data *data);
void write_fields(Sim_data *data, char *filename, char t);
void set_boundary(Sim_data *data);
void switch_n(Sim_data *data);
double divergence(Sim_data *data, int i, int j);
void update_pressure(Sim_data *data);
void mass_flow_condition(Sim_data *data);

#endif
