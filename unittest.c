#include "poisson.h"
#include "functions.h"

int main(int argc, char *argv[]) {
    int N = 2;
    int M = 2;

    double *u = calloc((N+2)*(M+1), sizeof(double));
    double *v = calloc((N+1)*(M+2), sizeof(double));
    double *u_pre = calloc((N+2)*(M+1), sizeof(double));
    double *v_pre = calloc((N+1)*(M+2), sizeof(double));
    double *P = calloc(N*M, sizeof(double));

    for (int i = 0; i < (N+2)*(M+1); i++) {
        u[i] = i;
        u_pre[i] = i;
    }
    for (int i = 0; i < (N+1)*(M+2); i++) {
        v[i] = i;
        v_pre[i] = i;
    }
    for (int i = 0; i < M*N; i++) {
        P[i] = i;
    }

    Sim_data *data = malloc(sizeof(Sim_data));
    data->u = u;
    data->v = v;
    data->u_pre = u_pre;
    data->v_pre = v_pre;
    data->P = P;
    data->M = M;
    data->N = N;
    data->nu = 1.0;
    data->h = 1.0;

    printf("%.1f should be = 25.0\n", get_a(u, v, M, N, 1.0, 1, 1));
    printf("%.1f should be = 14.5\n", get_b(u, v, M, N, 1.0, 1, 1));
    printf("%.1f should be = 2.0\n", get_Pu(P, M, N, 1.0, 1, 1));
    printf("%.1f should be = 1.0\n", get_Pv(P, M, N, 1.0, 1, 1));
    printf("%.1f should be = 0.0\n", get_lapu(u, M, N, 1.0, 1, 1));
    printf("%.1f should be = 0.0\n", get_lapv(v, M, N, 1.0, 1, 1));
    printf("RHS = %.1f should be = -27.0\n", get_RHSu(data, 1, 1));
    printf("RHS = %.1f should be = -15.5\n", get_RHSv(data, 1, 1));
}
