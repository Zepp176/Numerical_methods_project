#include <math.h>
#include "poisson.h"

double get_u(double *data, int M, int N, float idx_i, float idx_j) {
    int i = (int) idx_i;
    int j = (int) idx_j;
    char semi_i = !(idx_i == i);
    char semi_j = !(idx_j == j);

    int idx = 1 + (j-1) + (N+2)*(i-1);

    // Si au milieu de la cellule
    if (!semi_i && !semi_j) {
        return (data[idx]+data[idx+N+2])/2;
    }

    // Si en haut à droite de la cellule
    else if (semi_i && semi_j) {
        return (data[idx+N+2]+data[idx+N+3])/2;
    }

    // Si à droite
    else if (semi_i) {
        return data[idx+N+2];
    }

    return (data[idx]+data[idx+1]+data[idx+N+2]+data[idx+N+3])/4;
}

double get_v(double *data, int M, int N, float idx_i, float idx_j) {
    int i = (int) idx_i;
    int j = (int) idx_j;
    char semi_i = !(idx_i == i);
    char semi_j = !(idx_j == j);

    int idx = (j-1) + (N+1)*i;

    // Si au milieu de la cellule
    if (!semi_i && !semi_j) {
        return (data[idx]+data[idx+1])/2;
    }

    // Si en haut à droite de la cellule
    else if (semi_i && semi_j) {
        return (data[idx+1]+data[idx+N+2])/2;
    }

    // Si à droite
    else if (semi_i) {
        return (data[idx]+data[idx+1]+data[idx+N+1]+data[idx+N+2])/4;
    }

    return data[idx+1];
}

double get_a(double *u, double *v, int M, int N, double h, float idx_i, float idx_j) {
    return get_u(u, M, N, idx_i, idx_j)*(get_u(u, M, N, idx_i+1.0, idx_j) - get_u(u, M, N, idx_i-1.0, idx_j))/(2*h)
         + get_v(u, M, N, idx_i, idx_j)*(get_u(u, M, N, idx_i, idx_j+1.0) - get_u(u, M, N, idx_i, idx_j-1.0))/(2*h);
}

double get_b(double *u, double *v, int M, int N, double h, float idx_i, float idx_j) {
    return get_u(u, M, N, idx_i, idx_j)*(get_v(u, M, N, idx_i+1.0, idx_j) - get_v(u, M, N, idx_i-1.0, idx_j))/(2*h)
         + get_v(u, M, N, idx_i, idx_j)*(get_v(u, M, N, idx_i, idx_j+1.0) - get_v(u, M, N, idx_i, idx_j-1.0))/(2*h);
}

double get_grad_P(double *P, int M, int N, float idx_i, float idx_j) {
    int i = (int) idx_i;
    int j = (int) idx_j;
    char semi_i = !(idx_i == i);
    char semi_j = !(idx_j == j);

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
