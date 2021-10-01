#include "en-tête.h"
#ifndef fonction_
#define fonction_
#include "fonction_coeur.c"
#include "matrice_test.c"
#endif

// Jacobi

bool DiagDominante(double** MatriceA, int dim) {
    bool test = true;
    double somme1 = 0, somme2 = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            if (i == j)
                somme1 = MatriceA[i][j];
            else
                somme2 += MatriceA[i][j];
        }
        if (somme1 <= somme2 || somme1 == 0) test = false;
        somme1 = 0;
        somme2 = 0;
    }
    return test;
}

double** E(double** MatriceA, int dim) {
    double** E = (double**)malloc(dim * sizeof(double*));
    for (int i = 0; i < dim; i++) E[i] = (double*)malloc(dim * sizeof(int));

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            if (i > j) E[i][j] = -MatriceA[i][j];
    return E;
}

double** F(double** MatriceA, int dim) {
    double** F = (double**)malloc(dim * sizeof(double*));
    for (int i = 0; i < dim; i++) F[i] = (double*)malloc(dim * sizeof(int));

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            if (i < j) F[i][j] = -MatriceA[i][j];
    return F;
}

double** D1(double** MatriceA, int dim) {
    double** D1 = (double**)malloc(dim * sizeof(double*));
    for (int i = 0; i < dim; i++) D1[i] = (double*)malloc(dim * sizeof(int));

    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            if (i == j && MatriceA[i][j] != 0) D1[i][j] = 1 / (MatriceA[i][j]);
    return D1;
}

void calculD1EF(double** M, double** MatriceA, int dim) {
    // création variables
    double** e = E(MatriceA, dim);
    double** f = F(MatriceA, dim);
    double** d1 = D1(MatriceA, dim);
    // stocke la matrice après la somme
    double** S = creationMatriceVide(dim);
    // somme
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) S[i][j] = e[i][j] + f[i][j];
    // produit
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
            for (int a = 0; a < dim; a++) M[i][j] += d1[i][a] * S[a][j];
}

void calculD1b(double* M, double** MatriceA, double* solution, int dim) {
    // création variable
    double** d1 = D1(MatriceA, dim);

    // produit
    for (int i = 0; i < dim; i++)
        for (int a = 0; a < dim; a++) {
            M[i] += d1[i][a] * solution[a];
            // affiche(M, dim);
        }
}

double* Jacobi(double** MatriceA, int dim, double* solution) {
    double* D1b = (double*)malloc(dim * sizeof(double));
    double* inconnu = (double*)malloc(dim * sizeof(double));
    double* precedent = (double*)malloc(dim * sizeof(double));
    double** D1EF = (double**)malloc(dim * sizeof(double*));
    for (int c = 0; c < dim; c++)
        D1EF[c] = (double*)malloc(dim * sizeof(double));

    calculD1EF(D1EF, MatriceA, dim);
    // afficheM(D1EF, dim);
    calculD1b(D1b, MatriceA, solution, dim);
    // affiche(D1b, dim);
    double seuil = 0.000001;
    do {
        for (int i = 0; i < dim; i++) {
            precedent[i] = inconnu[i];
            inconnu[i] = 0.0;
        }
        for (int i = 0; i < dim; i++) {
            for (int a = 0; a < dim; a++)
                inconnu[i] += D1EF[i][a] * precedent[a];
            inconnu[i] += D1b[i];
        }
    } while (absdouble(inconnu[0] - precedent[0]) > seuil);

    return inconnu;
}
