#include "en-tête.h"
#ifndef fonction_
#define fonction_
#include "fonction_coeur.c"
#include "matrice_test.c"
#endif

// PIVOT DE GAUSS

double** matriceAug(double** matrice, int dim, double* solution) {
    // création amtrice augmentée
    double** matriceA = (double**)malloc(dim * sizeof(double*));
    for (int i = 0; i < dim; i++)
        matriceA[i] = (double*)malloc((dim + 1) * sizeof(double));
    // remplissage matrice
    for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) matriceA[i][j] = matrice[i][j];
    for (int i = 0; i < dim; i++) matriceA[i][dim] = solution[i];
    return matriceA;
}

void inverseLignes(double** Matrice, int l1, int l2, int dim) {
    int a;
    for (int i = 0; i < dim; i++) {
        a = Matrice[l1][i];
        Matrice[l1][i] = Matrice[l2][i];
        Matrice[l2][i] = a;
    }
}

void opLignes(double** Matrice, int l1, int l2, int dim, double coef) {
    for (int i = 0; i < dim + 1; i++) Matrice[l1][i] += coef * Matrice[l2][i];
}

void echelonnage(double** M, int dim, int i) {
    for (int j = i + 1; j < dim; j++) {
        double coef = -(M[j][i] / M[i][i]);
        opLignes(M, j, i, dim, coef);
    }
}

double* Gauss(double** MatriceA, int dim, bool* test) {
    //échelonnage matrice
    for (int i = 0; i < dim; i++) {
        if (MatriceA[i][i] != 0.0)
            echelonnage(MatriceA, dim, i);
        else {
            double a = 0;
            int tour = i;
            while (a == 0 && tour < dim) {
                a = MatriceA[tour][i];
                tour++;
            }
            if (tour < dim) {
                inverseLignes(MatriceA, i, tour, dim + 1);
                echelonnage(MatriceA, dim, i);
            } else {
                *test = false;
                return NULL;
            }
        }
    }

    // résolution système
    double* resultat = (double*)malloc(dim * sizeof(double));
    resultat[dim - 1] = MatriceA[dim - 1][dim] / MatriceA[dim - 1][dim - 1];
    for (int i = dim - 2; i >= 0; i--) {
        double somme = 0;
        for (int j = i + 1; j < dim; j++) somme += MatriceA[i][j] * resultat[j];
        resultat[i] = (1 / MatriceA[i][i]) * (MatriceA[i][dim] - somme);
    }
    return resultat;
}
