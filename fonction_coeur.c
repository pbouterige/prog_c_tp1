#include "en-tête.h"

double** creationMatriceVide(int dim) {
    double** M = (double**)malloc(dim * sizeof(double*));
    for (int i = 0; i < dim; i++) M[i] = (double*)malloc(dim * sizeof(int));
    return M;
}

void afficheM(double** M, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) printf("%f ", M[i][j]);
        puts("");
    }
    puts("");
}

void afficheMA(double** M, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim + 1; j++) printf("%f ", M[i][j]);
        puts("");
    }
    puts("");
}

void affiche(double* M, int h) {
    for (int i = 0; i < h; i++) printf("%f ", M[i]);
    puts("");
}

void remplirSol(double* M, int dim, int test) {
    double a = 0.0;
    if (test == 1)
        for (int i = 0; i < dim; i++) M[i] = 1;
    else
        for (int i = 0; i < dim; i++) {
            printf("valeur en position %d :", i);
            scanf("%lf", &a);
            M[i] = a;
        }
}

void remplirM(double** M, int dim) {
    double a = 1;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            printf("valeur en %d, %d :", i, j);
            scanf("%lf", &a);
            M[i][j] = a;
        }
    }
}

double** matriceCreuse(int dim) {
    double** MC = creationMatriceVide(dim);
    int a;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            a = rand() % 10;
            if (a <= 7)
                MC[i][j] = 0.0;
            else
                MC[i][j] = (double)(rand() % 4 + 1);
        }
    }
    return MC;
}

double absdouble(double i) {
    if (i > 0)
        return i;
    else
        return -i;
}

void free2D(double** m, int dim) {
    for (int i = 0; i < dim; i++) free(m[i]);
    free(m);
}

void libereMemoire(double** matrice, double* resultat, double* solution,
                   int dim) {
    free2D(matrice, dim);
    free(resultat);
    free(solution);
}