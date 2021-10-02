#include "Gauss.c"
#include "Jacobi.c"
#include "en-tête.h"

typedef enum { REMPLIR = 1, TEST = 2, CREUSE = 3 } Matrix_type;

int main() {
    srand(time(NULL));
    int dim = 500;

    double** matrice = Bord(dim);
    double* solution = (double*)malloc(dim * sizeof(double));
    remplirSol(solution, dim, 1);

    puts("");
    puts("Avec la matrice A :\n");
    afficheM(matrice, dim);

    double* resultatGauss = NULL;
    double* resultatJacobi = NULL;
    double** matriceA = matriceAug(matrice, dim, solution);
    bool test = true;
    resultatGauss = Gauss(matriceA, dim, &test);

    if (DiagDominante(matrice, dim)) {
        int a = 1;
        resultatJacobi = Jacobi(matrice, dim, solution);
        puts("On a les solutions suivantes avec Jacobi:\n ");
        for (int i = 0; i < dim; i++) {
            printf("x%d = %.3f  ", i + 1, resultatJacobi[i]);
            if (a % 5 == 0) puts("");
            a++;
        }
        puts("");

    } else {
        puts(
            "désolé, votre matrice n'est pas diagonale dominante donc "
            "Jacobi n'est pas adapté.");
    }

    puts("");
    if (test) {
        int a = 1;
        puts("on a les solutions suivantes avec Gauss:\n ");
        for (int i = 0; i < dim; i++) {
            printf("x%d = %.3f  ", i + 1, resultatGauss[i]);
            if (a % 5 == 0) puts("");
            a++;
        }
        puts("");
    } else
        puts("Nous n'avons pas de solution à vous proposer");

    free2D(matriceA, dim);
    free2D(matrice, dim);
    free(solution);
    if (resultatGauss) free(resultatGauss);
    if (resultatJacobi) free(resultatJacobi);

    return 0;
}

/*
int main() {
    srand(time(NULL));
    int dim;
    puts("Quelle est la dimension de votre matrice ?");
    scanf("%d", &dim);
LA:
    puts("");
    puts(
        "Voulez-vous : \n Remplir votre matrice (1)\n Utiliser une matrice "
        "test (2)\n Utiliser une matrice creuse (3)\n");
    int choix;
    scanf("%d", &choix);
    double** matrice = creationMatriceVide(dim);
    if (choix == 1)
        remplirM(matrice, dim);
    else if (choix == 2) {
        int choixmt = 0;
    ICI:
        puts(
            " Quelle matrice test voulez-vous utiliser :\n 1)Bord       "
            "2)Dingdong   "
            "3)Franc\n 4)Hilberb_m  5)Hilberd_p  6)kms\n 7)Lehmer 8)Lotkin "
            "    "
            "9)Moler\n");
        scanf("%d", &choixmt);
        switch (choixmt) {
            case 1:
                matrice = Bord(dim);
                break;
            case 2:
                matrice = DingDong(dim);
                break;
            case 3:
                matrice = Franc(dim);
                break;
            case 4:
                matrice = Hilbert_m(dim);
                break;
            case 5:
                matrice = Hilbert_p(dim);
                break;
            case 6:
                puts("valeur de p?");
                double p;
                scanf("%lf", &p);
                matrice = kms(dim, p);
                break;
            case 7:
                matrice = Lehmer(dim);
                break;
            case 8:
                matrice = Lotkin(dim);
                break;
            case 9:
                matrice = Moler(dim);
                break;
            default:
                puts("désolé le nombre que vous avez entrez n'est pas
valide."); goto ICI; break;
        }
    } else if (choix == 3)
        matrice = matriceCreuse(dim);
    else {
        puts("désolé le nombre que vous avez entrez n'est pas valide.");
        goto LA;
    }
    puts("");
    double* solution = (double*)malloc(dim * sizeof(double));
SOL:
    puts("quelle est votre matrice résultat ? \n");
    puts("1) (1 , 1,...,1)");
    puts("2) remplir");
    int r = 0;
    scanf("%d", &r);
    if (r == 1)
        remplirSol(solution, dim, 1);
    else if (r == 2)
        remplirSol(solution, dim, 0);
    else {
        puts("désolé le nombre que vous avez entrez n'est pas valide.");
        goto SOL;
    }

GJ:
    puts("");
    puts(
        "voulez-vous utiliser la méthode de Gauss(g) ou de Jacobi(j) pour "
        "resoudre le systeme ?\n");
    char a;
    getchar();
    scanf("%c", &a);
    double* resultat = (double*)malloc(dim * sizeof(double));
    bool test = true;
    if (a == 'g') {
        double** matriceA = matriceAug(matrice, dim, solution);
        resultat = Gauss(matriceA, dim, &test);
        // free2D(matriceA, dim);
    } else if (a == 'j') {
        if (DiagDominante(matrice, dim) == true)
            resultat = Jacobi(matrice, dim, solution);
        else {
            puts(
                "désolé, votre matrice n'est pas diagonale dominante donc "
                "Jacobi n'est pas adapté.");
            goto GJ;
        }
    } else {
        puts("entrez une valeur correcte s'il vous plaît.");
        goto GJ;
    }
    puts("");
    puts("Avec la matrice A :\n");
    afficheM(matrice, dim);
    puts("Et la matrice résultat b :\n");
    affiche(solution, dim);
    puts("");
    if (test == true) {
        int a = 1;
        puts("on a les solutions suivantes :\n ");
        for (int i = 0; i < dim; i++) {
            printf("x%d = %f  ", i + 1, resultat[i]);
            if (a % 3 == 0) puts("");
            a++;
        }
        puts("");
    } else
        puts("Nous n'avons pas de solution à vous proposer");
}
 */