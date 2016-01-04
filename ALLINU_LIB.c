#include "ALLINU_LIB.h"
#include <stdio.h>
#include "cblas.h"
#include <math.h>
#include <stdlib.h>
#include "f2c.h"

double* matPD(double* A, int m, int n, double* B, int q){
  int i, j, k;
  double *C = (double*) calloc(m*q,sizeof(double));
    for(i=0;i<m;i++) for(j=0;j<q;j++) for(k=0;k<n;k++) C[i+m*j]+=A[i+m*k]*B[k+n*j];
    return C;
}

double* gaxPD(double* A, int m, int n, double* b){
  int i, j;
  double *c = (double*) calloc(m,sizeof(double));
    for(i=0;i<m;i++) for(j=0;j<n;j++) c[i]+=A[i+m*j]*b[j];
    return c;
}

double* extraireQ(int m, int n, double *A, double *tau ){
  int i,j,k,l,z;
  double *Q, *H1,*H2, *v;

    Q  = (double *) calloc(m*m, sizeof(double));
    v  = (double *) calloc(m,   sizeof(double));
    H1 = (double *) calloc(m*m, sizeof(double));
    H2 = (double *) calloc(m*m, sizeof(double));

    for(i=0;i<m;i++) for(j=0;j<m;j++) Q[i+j*m]=0.0, H1[i+j*m]=0.0, H2[i+j*m]=0.0; 

    k = 0, v[0] = 1;
    for (i=1;i<m;i++) v[i] = A[i+m*k];
  
    for(i=0; i<m;i++) for(j=0;j<m;j++) H1[i+m*j] = (i==j)-tau[0]*v[i]*v[j];
 
    for(k=1;k<n;k++){
        for(i=0;i<k;i++) v[i] = 0.0;
        v[k] = 1;
        for(i=k+1;i<m;i++) v[i] = A[i+m*k];
        for(i=0;i<m;i++) for(j=0;j<m;j++) H2[i+m*j] = (i==j)- tau[k]*v[i]*v[j];
        for(i=0;i<m;i++) for(j=0;j<m;j++) for(l=0;l<m;l++) Q[i+j*m] += H1[i+l*m]*H2[l+j*m];
        for(i=0;i<m;i++) for(j=0;j<m;j++) H1[i+j*m] = Q[i+j*m],Q[i+j*m] = 0.0;	
    }

    free(Q); free(H2); free(v);

    return H1;
}

void  aasen(int n, double *A, double *T, double *L, double *P) {

/* Declaration des variables */

  long int i, j, p, k, l, ind, jpiv;
  double tau, tempv, temp;
  double  *v, *h,*tempL, *tempP, *tempA, *tempAbis;

 
 /* Allocation de memoire */
 
    tempL = (double *)   calloc(n, sizeof( double   ) );
    tempP = (double *)   calloc(n, sizeof( double   ) );
    tempA = (double *)   calloc(n, sizeof( double   ) );
    tempAbis = (double *)   calloc(n, sizeof( double   ) );
    h     = (double *)   calloc(n, sizeof( double   ) );
    v     = (double *)   calloc(n, sizeof( double   ) ); 
 

printf("Aasen strategie \n");

/*Initialisation : L, P matrice identite */

 for( i = 0; i < n; i++ ){
   L[i*n+i] = 1; 
   P[i*n+i] = 1; 
      }

 /* Boucle */

 for( j = 0; j < n; j++ ){
   if (j==0) {
     T[0]=A[0]; 
     printf( "v = \n" );
     for (i=1; i < n; i++ ){
       v[i]=A[i+0*n]; 
     }
   } //end for j==1

 if (j!=0) {

     h[0]= T[1+0*n] * L[j+1*n]; 

     for (k = 1; k < j; k++){
       h[k]=T[k+(k-1)*n] * L[j+(k-1)*n] + T[k+k*n]* L[j+k*n] + T[(k+1)+k*n]*L[j+(k+1)*n];
     } 
   

    h[j]=A[(j)+(j)*n];
     for (k = 0; k < j; k++){ 
      h[j]=h[j] - L[j+k*n] * h[k];
     } 

     T[j+j*n] = h[j] - T[j+(j-1)*n]* L[j+(j-1)*n];

     for (k=j+1; k < n; k++){
       temp=0;
     for (l=0; l < j+1; l++){  
     temp= temp+L[k+l*n]*h[l];
	 }

       v[k]=A[k+j*n]-temp;
     }

   } // end for j!=1


if (j<=n-2){

     ind=j;
     ind=j+1;
     tau=v[j+1];
   for (i=j+2; i < n; i++){ 
     if (abs(tau) < abs(v[i])) {ind=i; tau=v[i];}
    }
 
   jpiv=ind;
 
   //echange valeurs v
  
   tempv=v[jpiv];
   v[jpiv]=v[j+1];
   v[j+1]=tempv;
    
   
   //echange valeurs L
   for (l=1; l < j+1; l++){ 
     tempL[l]= L[(jpiv)+l*n];
     L[(jpiv)+l*n]=L[(j+1)+l*n];
         L[(j+1)+l*n]=tempL[l];
      }

   //echange valeurs P
   for (l=0; l < n; l++){ 
     tempP[l]= P[(jpiv)+l*n];
     P[(jpiv)+l*n]=P[(j+1)+l*n];
     P[(j+1)+l*n]=tempP[l];
      }

   //echange valeurs A
   for (l=j; l < n; l++){ 
     tempA[l]= A[(jpiv)+l*n];
     A[(jpiv)+l*n]=A[(j+1)+l*n];
     A[(j+1)+l*n]=tempA[l];
      }

     for (l=j; l < n; l++){ 
       tempAbis[l]= A[l+(jpiv)*n];
       A[l+(jpiv)*n]=A[l+(j+1)*n];
         A[l+(j+1)*n]=tempAbis[l];
      }

   T[j+1+j*n]=v[j+1];
   T[j+(j+1)*n]=v[j+1];
     
   }// end j<=n-2 

 if ((tau>0) && (j<=n-3)) {
   for (l=j+1; l < n; l++){ 
     L[l+(j+1)*n]= v[l]/v[j+1]; 
      }
   } //end if tau 


 }//end boucle


} //end Aasen

