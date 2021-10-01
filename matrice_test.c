#include "en-tête.h"

double** initM(int dim) {
    double** M = (double**)malloc(dim * sizeof(double*));
    for (int i = 0; i < dim; i++) M[i] = (double*)malloc(dim * sizeof(double));
    return M;
}

double** Bord(int dim) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (i == j)
                M[i][j] = 1;
            else if (i == 0 || j == 0) {
                if (i == 0)
                    M[i][j] = pow(2, -j);
                else
                    M[i][j] = pow(2, -i);
            } else
                M[i][j] = 0;
        }
    }
    return M;
}

double** DingDong(int dim) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) M[i][j] = 1 / (2 * (dim - i - j + 3.5));
    return M;
}

double** Franc(int dim) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            if (i >= j + 2)
                M[i][j] = 0;
            else
                M[i][j] = fmin(j + 1, i + 1);
        }
    return M;
}

double** Hilbert_m(int dim) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) M[i][j] = (1 / (double)(i + j + 1));
    return M;
}

double** Hilbert_p(int dim) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) M[i][j] = 1 / (double)(i + j + 3);
    return M;
}

double** kms(int dim, double p) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) M[i][j] = pow(p, abs(i - j));
    return M;
}

double** Lehmer(int dim) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            if (i <= j)
                M[i][j] = (double)(i + 1) / (double)(j + 1);
            else
                M[i][j] = (double)(j + 1) / (double)(i + 1);
        }
    return M;
}

double** Lotkin(int dim) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            if (i == 0)
                M[i][j] = 1;
            else
                M[i][j] = 1 / (double)(i + j + 1);
        }
    return M;
}

double** Moler(int dim) {
    double** M = initM(dim);
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
            if (i == j)
                M[i][j] = i + 1;
            else
                M[i][j] = fmin(i + 1, j + 1) - 2;
        }
    return M;
}