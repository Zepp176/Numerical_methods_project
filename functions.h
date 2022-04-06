#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

double get_u(double *data, int M, int N, float idx_i, float idx_j);
double get_v(double *data, int M, int N, float idx_i, float idx_j);
double get_a(double *u, double *v, int M, int N, double h, float idx_i, float idx_j);
double get_b(double *u, double *v, int M, int N, double h, float idx_i, float idx_j);
double get_grad_P(double *P, int M, int N, float idx_i, float idx_j);

#endif
