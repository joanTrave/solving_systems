
/* Autor: Big Vuk (@redKuv) para Informaticatech
 * 
 * Funciones de soporte para los apuntes de Algebra Lineal implementadas.
 * 
 * INDICE (por orden de implementación):
 * 
 * void resTinf (int n, double **L, double *x, double *b)
 * void resTsup (int n, double **L, double *x, double *b)
 * void prodMatVec (int m, int n, double **A, double *x, double *y)
 * void prodMatMat (int m, int n, int p, double **A, double **B, double **C)
 * 
 * Cualquier duda pueden contactar con el autor.
 * 
 * Este documento puede ser compartido libremente con fines educativos
 * o divulgativos.
 */


#include <stdio.h>
#include <stdlib.h>
#include "funciones_soporte_al.h"


void resTinf (int n, double **L, double *x, double *b){
	
	/*Resuelve el sistema Lx = b, donde:
	 * L es una matriz triangular inferior
	 * x el vector incognito
	 * b el vector independiente
	 * n las dimensiones de la matriz
	 */
	
    int i, j;	/* Variables auxiliares */
    x[0] = b[0];
    for (i = 1; i < n; i++){
        x[i] = b[i];
        for (j = 0; j < i; j++){
            x[i] -= L[i][j]*x[j];
        }
    }
}

void resTsup (int n, double **L, double *x, double *b){
	
	/* Resuelve el sistema Lx = b, donde:
	 * L es una matriz triangular superior
	 * x el vector incognito
	 * b el vector independiente
	 * n las dimensiones de la matriz
	 */
	
    int i, j;	/* Variables auxiliares */
    x[n-1] = b[n-1];
    for (i = n-2; i >= 0; i--){
        x[i] = b[i];
        for (j = i+1; j < n; j++){
            x[i] -= L[i][j]*x[j];
        }
    }
}

void prodMatVec (int m, int n, double **A, double *x, double *y){
	
	/* Realiza la multiplicacion entre una matriz de dimensiones m*n
	 * y un vector de dimension n:
	 * A*x = y
	 */ 
	
    int i, j;	/* Variables auxiliares */
    for (i = 0; i < m; i++){
        y[i] = 0;
        for (j = 0; j < n; j++){
            y[i] += A[i][j]*x[j];
        }
    }
}

void prodMatMat (int m, int n, int p, double **A, double **B, double **C){
	
	/* Realiza la multiplicacion entre dos matrices A y B
	 * A <- m*n     B <- n*m
	 * C = A*B <- m*m
	 */
	
    int i, j, k;	/* Variables auxiliares */
    for (i = 0; i < m; i++){
        for (j = 0; j < p; j++){
            C[i][j] = 0;
            for (k = 0; k < n; k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}
