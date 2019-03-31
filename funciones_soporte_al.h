
/* Autor: Big Vuk (@redKuv) para Informaticatech
 * 
 * Funciones de soporte para los apuntes de Algebra Lineal declaradas.
 * 
 * Cualquier duda pueden contactar con el autor.
 * 
 * Este documento puede ser compartido libremente con fines educativos
 * o divulgativos.
 */

#include <stdio.h>
#include <stdlib.h>


void resTinf (int n, double **L, double *x, double *b);
void resTsup (int n, double **L, double *x, double *b);
void prodMatVec (int m, int n, double **A, double *x, double *y);
void prodMatMat (int m, int n, int p, double **A, double **B, double **C);
