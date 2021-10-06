#include "Gauss.c"
#include "Jacobi.c"
#include "en-tête.h"

int main() {
    srand(time(NULL));
    int dim = 500;

    double** matrice = matriceCreuse(dim);
    double* solution = (double*)malloc(dim * sizeof(double));
    remplirSol(solution, dim, 1);

    // puts("");
    // puts("Avec la matrice A :\n");
    // afficheM(matrice, dim);

    double* resultatGauss = NULL;
    double* resultatJacobi = NULL;
    double* verif_Gauss = NULL;
    double* verif_Jacobi = NULL;
    double** matriceA = matriceAug(matrice, dim, solution);
    bool test = true;
    int temps_debut, temps_fin;
    if (DiagDominante(matrice, dim)) {
        int itération = 0;
        temps_debut = clock();
        resultatJacobi = Jacobi(matrice, dim, solution, &itération);
        temps_fin = clock();
        puts("\n                 JACOBI                ");
        verif_Jacobi = multMatrice(dim, matrice, resultatJacobi);
        puts("-----------------------------------------");
        printf("taille matrice     : %d x %d\n", dim, dim);
        pourcentage_ecart(verif_Jacobi, dim);
        printf("nombres d'itérations : %d", itération);
        printf("\ntemps execution     : %f s\n",
               (float)(temps_fin - temps_debut) / 1000000);
        fonction_erreur(verif_Jacobi, dim);
        puts("-----------------------------------------");
    } else {
        puts(
            "désolé, votre matrice n'est pas diagonale dominante donc "
            "Jacobi n'est pas adapté.");
    }
    temps_debut = clock();
    resultatGauss = Gauss(matriceA, dim, &test);
    temps_fin = clock();
    if (test) {
        puts("\n                 GAUSS                ");
        verif_Gauss = multMatrice(dim, matrice, resultatGauss);
        puts("-----------------------------------------");
        printf("taille matrice     : %d x %d\n", dim, dim);
        pourcentage_ecart(verif_Gauss, dim);
        printf("temps execution        : %f s\n",
               (float)(temps_fin - temps_debut) / 1000000);
        fonction_erreur(verif_Gauss, dim);
        puts("-----------------------------------------");
    } else
        puts("Nous n'avons pas de solution à vous proposer");

    free2D(matriceA, dim);
    free2D(matrice, dim);
    free(solution);
    if (resultatGauss) free(resultatGauss);
    if (resultatJacobi) free(resultatJacobi);
    if (verif_Jacobi) free(verif_Jacobi);
    if (verif_Gauss) free(verif_Gauss);

    return 0;
}
