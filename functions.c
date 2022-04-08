#include <math.h>
#include "poisson.h"
#include "functions.h"

double get_u(Sim_data *data, int t, float idx_i, float idx_j) {
    int N = data->N;
    double *u;
    if (t == -1) {u = data->u_pre;}
    else if (t == 1) {u = data->u_star;}
    else {u = data->u;}

    int i = (int) idx_i;
    int j = (int) idx_j;
    char semi_i = !(idx_i == i);
    char semi_j = !(idx_j == j);

    int idx = 1 + (j-1) + (N+2)*(i-1);

    // Si au milieu de la cellule
    if (!semi_i && !semi_j) {
        return (u[idx]+u[idx+N+2])/2;
    }

    // Si en haut à droite de la cellule
    else if (semi_i && semi_j) {
        return (u[idx+N+2]+u[idx+N+3])/2;
    }

    // Si à droite
    else if (semi_i) {
        return u[idx+N+2];
    }

    return (u[idx] + u[idx+1] + u[idx+N+2] + u[idx+N+3]) / 4;
}

double get_v(Sim_data *data, int t, float idx_i, float idx_j) {
    int N = data->N;
    double *v;
    if (t == -1) {v = data->v_pre;}
    else if (t == 1) {v = data->v_star;}
    else {v = data->v;}

    int i = (int) idx_i;
    int j = (int) idx_j;
    char semi_i = !(idx_i == i);
    char semi_j = !(idx_j == j);

    int idx = (j-1) + (N+1)*i;

    // Si au milieu de la cellule
    if (!semi_i && !semi_j) {
        return (v[idx]+v[idx+1])/2;
    }

    // Si en haut à droite de la cellule
    else if (semi_i && semi_j) {
        return (v[idx+1]+v[idx+N+2])/2;
    }

    // Si à droite
    else if (semi_i) {
        return (v[idx] + v[idx+1] + v[idx+N+1] + v[idx+N+2]) / 4;
    }

    return v[idx+1];
}

double get_a(Sim_data *data, int t, float i, float j) {
    double h = data->h;

    return get_u(data, t, i, j)*(get_u(data, t, i+1.0, j) - get_u(data, t, i-1.0, j))/(2*h)
         + get_v(data, t, i, j)*(get_u(data, t, i, j+1.0) - get_u(data, t, i, j-1.0))/(2*h);
}

double get_b(Sim_data *data, int t, float i, float j) {
    double h = data->h;

    return get_u(data, t, i, j)*(get_v(data, t, i+1.0, j) - get_v(data, t, i-1.0, j))/(2*h)
         + get_v(data, t, i, j)*(get_v(data, t, i, j+1.0) - get_v(data, t, i, j-1.0))/(2*h);
}

double get_grad_P(Sim_data *data, float idx_i, float idx_j) {
    double h = data->h;
    double *P = data->P;
    int N = data->N;

    int i = (int) idx_i;
    int j = (int) idx_j;
    char semi_i = !(idx_i == i);

    int idx = j + i*N;

    // pour vitesse u
    if (semi_i) {
        return (P[idx+N] - P[idx]) / h;
    }

    // pour vitesse v
    else {
        return (P[idx+1] - P[idx]) / h;
    }
}

double get_laplacian(Sim_data *data, float idx_i, float idx_j) {
    double h = data->h;

    int i = (int) idx_i;
    char semi_i = !(idx_i == i);

    // pour vitesse u
    if (semi_i) {
        return (-4*get_u(data, 0, idx_i, idx_j) + get_u(data, 0, idx_i+1.0, idx_j) + get_u(data, 0, idx_i-1.0, idx_j)
                                                + get_u(data, 0, idx_i, idx_j+1.0) + get_u(data, 0, idx_i, idx_j-1.0)) / (h*h);
    }

    // pour vitesse v
    else {
        return (-4*get_v(data, 0, idx_i, idx_j) + get_v(data, 0, idx_i+1.0, idx_j) + get_v(data, 0, idx_i-1.0, idx_j)
                                                + get_v(data, 0, idx_i, idx_j+1.0) + get_v(data, 0, idx_i, idx_j-1.0)) / (h*h);
    }
}

void init_sim_data(Sim_data *data, int res, double Re) {
    int N = 5*res;
    int M = 15*res;

    data->N = N;
    data->M = M;

    data->nu = 1.0;
    data->H_box = 1.0;
    data->U_inf = Re;
    data->h = data->H_box/res;
    data->dt = 1.0/(4*res*res);

    data->u      = calloc((N+2)*(M+1), sizeof(double)); // u_n+1
    data->v      = calloc((N+1)*(M+2), sizeof(double));
    data->u_star = calloc((N+2)*(M+1), sizeof(double)); // u*
    data->v_star = calloc((N+1)*(M+2), sizeof(double));
    data->u_pre  = calloc((N+2)*(M+1), sizeof(double)); // u_n
    data->v_pre  = calloc((N+1)*(M+2), sizeof(double));
    data->P      = calloc(N*M, sizeof(double));
}

void free_sim_data(Sim_data *data) {
    free(data->u);
    free(data->v);
    free(data->u_pre);
    free(data->v_pre);
    free(data->u_star);
    free(data->v_star);
    free(data->P);
}

void write_fields(Sim_data *data, char *filename) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("ERROR opening %s\n", filename);
        exit(1);
    }

    int M = data->M;
    int N = data->N;
    double *u = data->u_star;
    double *v = data->v_star;
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
    }

    // Lateral boundaries
    for (int i = 1; i < M+1; i++) {
        u_star[(N+2)*i]       = u[(N+2)*i + 1]; // ghost points : no vorticity
        u_star[(N+2)*i + N+1] = u[(N+2)*i + N];
        v_star[(N+1)*i]     = 0.0; // no tangential velocity
        v_star[(N+1)*i + N] = 0.0;
    }

    // Right boundary
    int idx;
    for (int j = 1; j < N+1; j++) {
        idx = M*(N+2) + j;
        u_star[idx] = u[idx] - U_inf*(u[idx] - u[idx-(N+2)])/h;
    }

    double vort_p; double vort_m;
    for (int j = 1; j < N; j++) {
        idx = (M+1)*(N+1) + j;
        vort_p = get_v(data, 0, M+1, j+0.5) - get_v(data, 0, M, j+0.5) + get_u(data, 0, M+0.5, j) - get_u(data, 0, M+0.5, j+1);
        vort_m = get_v(data, 0, M, j+0.5) - get_v(data, 0, M-1, j+0.5) + get_u(data, 0, M-0.5, j) - get_u(data, 0, M-0.5, j+1);
        v_star[idx] = (1 - dt/h)*vort_p + dt/h * vort_m - get_u(data, 1, M+0.5, j) + get_u(data, 1, M+0.5, j+1) + get_v(data, 1, M, j+0.5);
    }
}
