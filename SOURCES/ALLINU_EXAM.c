#include "ALLINU_EXAM.h"

int main(){
  double *A, *A2, *b, *b2, *L, *T, *P, *D, *SD1, *SD2, *PT, *WORK;
  int SIZE, m, POSX, POSY, NUML, NMBR, i, j, NRHS=1, INFO, *IPIV, NUMM;
  char NAME[20], FNAM[30], UPLO='L', TRANS='N', DIAG='U';
  FILE *FILI,*FILI2,*FILO, *FIL3;

    FIL3 = fopen("INPUT.dat","r");
// LECTURE DU NOMBRE DE MATRICE A TRAITER
    fscanf(FIL3,"%d",&NMBR);
    for(NUMM = 0; NUMM<NMBR;NUMM++){
// LECTURE DU NOM DE LA MATRICE ACTUELLE ET OUVERTURE DES FICHIERS
        fscanf(FIL3,"%s",NAME);
        FILI = fopen(NAME,"r");
        if(FILI==NULL){
            printf("Fichier %s introuvable !!!\n",NAME);
            continue;
        }
        sprintf(HEAD,"rhs_%s.dat",NAME);
        FILI2 = fopen(FNAM,"r");
        sprintf(HEAD,"result_%s.dat",NAME);
        FILO = fopen(FNAM,"w+");

// ALLOCATION DES DIFFERENTS TABLEAUX
        fscanf(FILI,"%d %d %d",&SIZE, &m, &NUML);
        A    = (double *) calloc(SIZE*SIZE,sizeof(double));
        A2   = (double *) calloc(SIZE*SIZE,sizeof(double));
        b    = (double *) calloc(SIZE,sizeof(double));
		b2   = (double *) calloc(SIZE,sizeof(double));
        P    = (double *) calloc(SIZE*SIZE,sizeof(double));
        PT   = (double *) calloc(SIZE*SIZE,sizeof(double));    
        L    = (double *) calloc(SIZE*SIZE,sizeof(double));
        T    = (double *) calloc(SIZE*SIZE,sizeof(double));
        D    = (double *) calloc(SIZE,sizeof(double));
        SD1  = (double *) calloc(SIZE-1,sizeof(double));
        SD2  = (double *) calloc(SIZE-1,sizeof(double));
        WORK = (double *) calloc(SIZE,sizeof(double));
        IPIV = (int*)     calloc(SIZE,sizeof(int));    

// LECTURE DES VALEURS DE A
        for(i=0;i<NUML;i++){
            fscanf(FILI,"%d %d",&POSX, &POSY);
            fscanf(FILI,"%lf",A+POSX-1+SIZE*(POSY-1));
            A2[(POSX-1)+SIZE*(POSY-1)] = A[(POSX-1)+SIZE*(POSY-1)];
        }

// LECTURE DES VALEURS DE b
        for(i=0;i<SIZE;i++){
            fscanf(FILI2,"%lf",b+i);
            b2[i] = b[i];
        }

// METHODE DE AASEN
//   CALCUL DES MATRICES T,L ET P
        aasen(SIZE,A,T,L,P);

//   RESOLUTION DES AUTRES EQUATIONS
        b = gaxPD(P,SIZE,SIZE,b);
        for(i=0;i<SIZE;i++) D[i] = T[i+SIZE*i];
        for(i=0;i<SIZE-1;i++) SD1[i] = T[i+SIZE*(i+1)], SD2[i] = T[i+1+SIZE*i];
        for(i=0;i<SIZE;i++) for(j=0;j<SIZE;j++) PT[i+SIZE*j]=P[j+SIZE*i];

        dtrtrs(&UPLO, &TRANS, &DIAG, &SIZE, &NRHS, L, &SIZE, b, &SIZE, &INFO);
        dgtsv(&SIZE, &NRHS, SD1, D, SD2, b, &SIZE, &INFO);
        TRANS = 'T';
        dtrtrs(&UPLO,&TRANS,&DIAG,&SIZE,&NRHS,L,&SIZE,b,&SIZE,&INFO);

//   CALCUL DE LA SOLUTION FINALE
        b = gaxPD(PT,SIZE,SIZE,b);

// METHODE DE BUNCH-PARLETT
        dsytrf_(&UPLO, &SIZE, A2, &SIZE, IPIV, WORK, &SIZE, &INFO);
        dsytrs_(&UPLO, &SIZE, &NRHS, A2, &SIZE, IPIV, b2, &SIZE, &INFO);

// ECRITURE DANS LE FICHIER DE SORTIE
        for(i=0;i<SIZE;i++) fprintf(FILO,"%lf\t%lf\SIZE",b[i],b2[i]);

// LIBERATION DE LA MEMOIRE ET FERMETURE DES FICHIERS
        free(A);free(b);free(A2);free(b2);free(P);free(PT);free(L);free(T);free(D);free(SD1);
        free(SD2);
        fclose(FILI); fclose(FILI2); fclose(FILO);
    }
    return(0);
}
