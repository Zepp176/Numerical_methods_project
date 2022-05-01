#include <math.h>
#include "poisson.h"
#include "functions.h"

double get_a(double *u, double *v, int M, int N, double h, int i, int j) {
    double u_ip32   = u[(i+1)*(N+2) + j];
    double u_ip12   = u[i*(N+2) + j];
    double u_im12   = u[(i-1)*(N+2) + j];
    double u_jp1    = u[i*(N+2) + j+1];
    double u_jm1    = u[i*(N+2) + j-1];
    double v_jp12   = v[i*(N+1) + j];
    double v_jm12   = v[i*(N+1) + j-1];
    double v_i1jp12 = v[(i+1)*(N+1) + j];
    double v_i1jm12 = v[(i+1)*(N+1) + j-1];

    double u_r = (u_ip32 + u_ip12)/2.0;
    double u_l = (u_ip12 + u_im12)/2.0;
    double v_u = (v_jp12 + v_i1jp12)/2.0;
    double v_d = (v_jm12 + v_i1jm12)/2.0;

    double a = 1.0/2.0 * (u_r * (u_ip32 - u_ip12) + u_l * (u_ip12 - u_im12));
          a += 1.0/2.0 * (v_u * (u_jp1  - u_ip12) + v_d * (u_ip12 - u_jm1 ));

    return a/h;
}

double get_b(double *u, double *v, int M, int N, double h, int i, int j) {
    double v_l = v[(i-1)*(N+1) + j];
    double v_m = v[i*(N+1) + j];
    double v_r = v[(i+1)*(N+1) + j];
    double v_u = v[i*(N+1) + j+1];
    double v_d = v[i*(N+1) + j-1];
    double u_ul = u[(i-1)*(N+2) + j+1];
    double u_ur = u[i*(N+2) + j+1];
    double u_dl = u[(i-1)*(N+2) + j];
    double u_dr = u[i*(N+2) + j];

    double u_rm = (u_ur + u_dr)/2.0;
    double u_lm = (u_ul + u_dl)/2.0;
    double v_um = (v_u + v_m)/2.0;
    double v_dm = (v_m + v_d)/2.0;

    double b = 1.0/2.0 * (u_rm * (v_r - v_m) + u_lm * (v_m - v_l));
          b += 1.0/2.0 * (v_um * (v_u - v_m) + v_dm * (v_m - v_d));

    return b/h;
}

double get_gradu(double *P, int M, int N, double h, int i, int j) {
    return (P[i*N + j-1] - P[(i-1)*N + j-1]) / h;
}

double get_gradv(double *P, int M, int N, double h, int i, int j) {
    return (P[(i-1)*N + j] - P[(i-1)*N + j-1]) / h;
}

double get_lapu(double *u, int M, int N, double h, int i, int j) {
    double u_r = u[(i+1)*(N+2) + j];
    double u_m = u[i*(N+2) + j];
    double u_l = u[(i-1)*(N+2) + j];
    double u_u = u[i*(N+2) + j+1];
    double u_d = u[i*(N+2) + j-1];

    return (u_r + u_l + u_u + u_d - 4*u_m) / (h*h);
}

double get_lapv(double *v, int M, int N, double h, int i, int j) {
    double v_l = v[(i-1)*(N+1) + j];
    double v_m = v[i*(N+1) + j];
    double v_r = v[(i+1)*(N+1) + j];
    double v_u = v[i*(N+1) + j+1];
    double v_d = v[i*(N+1) + j-1];

    return (v_l + v_r + v_u + v_d - 4*v_m) / (h*h);
}

double get_RHSu(Sim_data *data, int i, int j) {
    int M = data->M;
    int N = data->N;
    double h = data->h;
    double nu = data->nu;

    double *u = data->u;
    double *v = data->v;
    double *u_pre = data->u_pre;
    double *v_pre = data->v_pre;
    double *P = data->P;

    double Hn  = get_a(u,     v,     M, N, h, i, j);
    double Hnm = get_a(u_pre, v_pre, M, N, h, i, j);
    double gradP = get_gradu(P, M, N, h, i, j);
    double lapl = get_lapu(u, M, N, h, i, j);

    return -1.0/2.0 * (3.0 * Hn - Hnm) - gradP + nu * lapl;
}

double get_RHSv(Sim_data *data, int i, int j) {
    int M = data->M;
    int N = data->N;
    double h = data->h;
    double nu = data->nu;

    double *u = data->u;
    double *v = data->v;
    double *u_pre = data->u_pre;
    double *v_pre = data->v_pre;
    double *P = data->P;

    double Hn  = get_b(u,     v,     M, N, h, i, j);
    double Hnm = get_b(u_pre, v_pre, M, N, h, i, j);
    double gradP = get_gradv(P, M, N, h, i, j);
    double lapl = get_lapv(v, M, N, h, i, j);

    return -1.0/2.0 * (3.0 * Hn - Hnm) - gradP + nu * lapl;
}

void compute_star(Sim_data *data) {
    int M = data->M;
    int N = data->N;
    double dt = data->dt;

    double *u_star = data->u_star;
    double *v_star = data->v_star;
    double *u = data->u;
    double *v = data->v;

    for (int i = 1; i < M; i++) { // Pour les u
        for (int j = 1; j < N+1; j++) {
            u_star[i*(N+2) + j] = u[i*(N+2) + j] + dt * get_RHSu(data, i, j);
        }
    }

    for (int i = 1; i < M+1; i++) { // Pour les v
        for (int j = 1; j < N; j++) {
            v_star[i*(N+1) + j] = v[i*(N+1) + j] + dt * get_RHSv(data, i, j);
        }
    }
}

void compute_next(Sim_data *data) {
    int M = data->M;
    int N = data->N;
    double dt = data->dt;
    double h = data->h;

    double *u = data->u;
    double *v = data->v;
    double *u_star = data->u_star;
    double *v_star = data->v_star;

    for (int i = 1; i < M; i++) { // Pour les u
        for (int j = 1; j < N+1; j++) {
            u[i*(N+2) + j] = u_star[i*(N+2) + j] - dt * get_gradu(data->phi, M, N, h, i, j);
        }
    }

    for (int i = 1; i < M+1; i++) { // Pour les v
        for (int j = 1; j < N; j++) {
            v[i*(N+1) + j] = v_star[i*(N+1) + j] - dt * get_gradv(data->phi, M, N, h, i, j);
        }
    }

    // Laterals of the box
    int res = data->res;
    for (int i = 1 + 3*res; i < 1 + 8*res; i++) {
        int j = 2*res;
        v[i*(N+1) + j]     = 0.0; // down
        v[i*(N+1) + j+res] = 0.0; // up
    }
    for (int j = 1 + 2*res; j < 1 + 3*res; j++) {
        int i = 3*res;
        u[i*(N+2) + j]         = 0.0; // left
        u[(i+5*res)*(N+2) + j] = 0.0; // right
    }
}

void init_sim_data(Sim_data *data, int res, double Re) {
    int N = 5*res;
    int M = 15*res;

    data->N = N;
    data->M = M;

    data->nu = 0.000001;
    data->H_box = 0.01;
    data->U_inf = Re*data->nu/data->H_box;
    data->h = data->H_box/res;
    data->dt = 1.0/(4.0*res*res);
    data->res = res;

    data->u      = calloc((N+2)*(M+1), sizeof(double)); // u_n+1
    data->v      = calloc((N+1)*(M+2), sizeof(double));
    data->u_star = calloc((N+2)*(M+1), sizeof(double)); // u*
    data->v_star = calloc((N+1)*(M+2), sizeof(double));
    data->u_pre  = calloc((N+2)*(M+1), sizeof(double)); // u_n
    data->v_pre  = calloc((N+1)*(M+2), sizeof(double));
    data->P      = calloc(N*M, sizeof(double));
    data->phi    = calloc(N*M, sizeof(double));
}

void free_sim_data(Sim_data *data) {
    free(data->u);
    free(data->v);
    free(data->u_pre);
    free(data->v_pre);
    free(data->u_star);
    free(data->v_star);
    free(data->P);
    free(data->phi);
}

void write_fields(Sim_data *data, char *filename, char t) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("ERROR opening %s\n", filename);
        exit(1);
    }

    int M = data->M;
    int N = data->N;
    double *u;
    double *v;
    if (t == 1) {
        u = data->u_star;
        v = data->v_star;
    } else if (t == -1) {
        u = data->u_pre;
        v = data->v_pre;
    } else {
        u = data->u;
        v = data->v;
    }
    double *P = data->P;

    fprintf(f, "%d %d\n", M, N);
    for (int i = 0; i < (N+2)*(M+1); i++) {
        fprintf(f, "%f ", u[i]);
    }
    fprintf(f, "\n");
    for (int i = 0; i < (N+1)*(M+2); i++) {
        fprintf(f, "%f ", v[i]);
    }
    fprintf(f, "\n");
    for (int i = 0; i < N*M; i++) {
        fprintf(f, "%f ", P[i]);
    }
    fclose(f);
}

void set_boundary(Sim_data *data) {
    int M = data->M;
    int N = data->N;
    int res = data->res;
    double U_inf = data->U_inf;
    double h = data->h;
    double dt = data->dt;
    double *u = data->u;
    double *v = data->v;
    double *u_star = data->u_star;
    double *v_star = data->v_star;

    // Left boundary
    for (int j = 1; j < N+1; j++) {
        u_star[j] = U_inf; // u = U_infinity
        v_star[j] = -1.0/5.0 * (v[j + 3*(N+1)] - 5*v[j + 2*(N+1)] + 15*v[j + 1*(N+1)]); // ghost points : no tangential velocity

        u[j] = u_star[j];
        v[j] = v_star[j];
    }

    // Lateral boundaries
    for (int i = 1; i < M; i++) {
        u_star[(N+2)*i]       = u[(N+2)*i + 1]; // ghost points : no vorticity
        u_star[(N+2)*i + N+1] = u[(N+2)*i + N];

        u[(N+2)*i] = u_star[(N+2)*i];
        u[(N+2)*i + N+1] = u_star[(N+2)*i + N+1];
    }
    for (int i = 1; i < M+1; i++) {
        v_star[(N+1)*i]     = 0.0; // no tangential velocity
        v_star[(N+1)*i + N] = 0.0;

        v[(N+1)*i] = v_star[(N+1)*i];
        v[(N+1)*i + N] = v_star[(N+1)*i + N];
    }

    // Right boundary
    int idx;
    for (int j = 1; j < N+1; j++) {
        idx = M*(N+2) + j;
        u_star[idx] = u[idx] - (u[idx] - u[idx-(N+2)]) * U_inf * dt / h;

        u[idx] = u_star[idx];
    }

    double vort_p; double vort_m; double RHS;
    for (int j = 1; j < N; j++) {
        idx = M*(N+2) + j;
        vort_p = u[idx] - u[idx+1];
        idx = (M-1)*(N+2) + j;
        vort_m = u[idx] - u[idx+1];
        idx = M*(N+1) + j;
        vort_p += v[idx+N+1] - v[idx];
        idx = (M-1)*(N+1) + j;
        vort_m += v[idx+N+1] - v[idx];

        RHS = vort_p - U_inf*dt/h * (vort_p - vort_m);
        v_star[(M+1)*(N+1) + j] = v_star[M*(N+1) + j] + u_star[M*(N+2) + j+1] - u_star[M*(N+2) + j] + RHS;

        v[(M+1)*(N+1) + j] = v_star[(M+1)*(N+1) + j];
    }

    // Laterals of the box
    for (int i = 1 + 3*res; i < 1 + 8*res; i++) {
        int j = 2*res;
        v_star[i*(N+1) + j]     = 0.0; // down
        v_star[i*(N+1) + j+res] = 0.0; // up
    }
    for (int i = 3*res + 1; i < 8*res; i++) {
        int j = 3*res;

        int idx = i*(N+2) + j;
        u_star[idx] = -1.0/5.0 * (15*u[idx+1] - 5*u[idx+2] + u[idx+3]);

        idx = i*(N+2) + j-res+1;
        u_star[idx] = -1.0/5.0 * (15*u[idx-1] - 5*u[idx-2] + u[idx-3]);
    }
    for (int j = 1 + 2*res; j < 1 + 3*res; j++) {
        int i = 3*res;
        u_star[i*(N+2) + j]         = 0.0; // left
        u_star[(i+5*res)*(N+2) + j] = 0.0; // right
    }
    for (int j = 2*res + 1; j < 3*res; j++) {
        int i = 3*res + 1;
        int idx = i*(N+1) + j;
        v_star[idx] = -1.0/5.0 * (15*v[idx-(N+1)] - 5*v[idx-2*(N+1)] + v[idx-3*(N+1)]);

        i = 8*res;
        idx = i*(N+1) + j;
        v_star[idx] = -1.0/5.0 * (15*v[idx+(N+1)] - 5*v[idx+2*(N+1)] + v[idx+3*(N+1)]);
    }
}

void switch_n(Sim_data *data) {
    int M = data->M;
    int N = data->N;
    double *u = data->u;
    double *v = data->v;
    double *u_pre = data->u_pre;
    double *v_pre = data->v_pre;

    for (int i = 0; i < (N+2)*(M+1); i++) {
        u_pre[i] = u[i];
    }
    for (int i = 0; i < (N+1)*(M+2); i++) {
        v_pre[i] = v[i];
    }
}

double divergence(Sim_data *data, int i, int j) {
    int N = data->N;
    double *u = data->u_star;
    double *v = data->v_star;
    return u[i*(N+2) + j] - u[(i-1)*(N+2) + j] + v[i*(N+1) + j] - v[i*(N+1) +j-1]; // ne pas mettre 1/h ici!!! déjà pris en compte dans A
}

void update_pressure(Sim_data *data) {
    int M = data->M;
    int N = data->N;
    double *P = data->P;
    double *phi = data->phi;

    for (int i = 0; i < M*N; i++) {
        P[i] += phi[i];
    }
}

void mass_flow_condition(Sim_data *data) {
    int M = data->M;
    int N = data->N;

    double sum_i = 0.0;
    double sum_o = 0.0;

    for (int j = 1; j < N+1; j++) {
        sum_i += data->u_star[j];
        sum_o += data->u_star[M*(N+2) + j];
    }

    for (int j = 1; j < N+1; j++) {
        data->u_star[M*(N+2) + j] += (sum_i - sum_o)/N;
    }
}
