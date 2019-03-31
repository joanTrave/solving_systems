
/* Autor: Joan Travé (@dzerkuv)
 * 
 * ----------------------------------------------------------------------------
 * 
 * Se resolverá aqui el sistema Ax = b, donde:
 * 
 * A es una matriz m*n, con m > n (m las ecuaciones y n las incognitas)
 * x es el vector con incognitas
 * b el vector independiente
 * 
 * ----------------------------------------------------------------------------
 * 
 * EXPLICACIÓN: Si tenemos planteado un sistema (Ax = b) con más ecuaciones que 
 * incognitas es muy fácil que sea éste un sistema indeterminado (queda a cargo 
 * del alumno entender eso). Así pues tendremos que aproximar la solucion, que 
 * podemos llamar x*. Donde la norma de Ax - b proxima a zero. Para eso 
 * transfomaremos la matriz A en  cuadrangular (de dimension n*n). Podemos ver 
 * que ésto ultimo lo conseguiremos multiplicando el sistema por la transpuesta 
 * de A (At).
 * 
 * Una vez tengamos el sistema AtAx = Atb, siendo AtA = At*A y Atb = At*b,
 * realizaremos una transformacion ldlt en la matriz AtA. Luego resolveremos el
 * sistema de la siguiente manera:
 * 
 * AtAx = (LDLt)*x = Atb --> L*z = Atb (encontramos z) --> D*y = z (encontramos
 * y), y, por ultimo Lt*x = y (de solucion trivial).
 * 
 * ---------------------------------------------------------------------------
 * 
 * INDICE (por orden de implementación):
 * 
 * En funciones_soporte_al.h:
 * 
 * void resTinf (int n, double **L, double *x, double *b)
 * void resTsup (int n, double **L, double *x, double *b)
 * void prodMatVec (int m, int n, double **A, double *x, double *y)
 * void prodMatMat (int m, int n, int p, double **A, double **B, double **C)
 * 
 * En este documento:
 * 
 * int ldlt (int n, double **A, double tol)
 * int main()
 * 
 * ----------------------------------------------------------------------------
 * 
 * Cualquier duda pueden contactar con el autor.
 * 
 * Este documento puede ser compartido libremente con fines educativos
 * o divulgativos.
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "funciones_soporte_al.h"

int ldlt (int n, double **A, double tol);

int main(){
	
	/* DECLARACION DE VARIABLES */
	int i, j;
	/* Variables auxiliares */
	int m, n;
	/* Dimensiones */
	double tolerancia, norma = 0;
	double *b, *Atb, *x, *y, *z, *Ax;
	/* b: vector independiente
	 * Atb: matriz A traspuesta multiplicada por el vector independiente
	 * x: respuesta
	 * y,z: vectores auxiliares para resolver una vez descompuesta
	 * Ax: matriz A * vector x, para calcular posteriormente la norma */
	double **A, **At, **AtA;			
	/* A:matriz con los coeficientes de las ecuaciones
	 * At: traspuesta de A
	 * AtA: Matriz trasposada d'A * matriz A */
	
	printf ("***************RESOLVER SISTEMAS*****************\n\n");
	
	/* Leemos las dimensiones de A y comprovamos que son realistas */
	printf ("Dimensiones de A: ");
	scanf ("%d %d", &m, &n);
	while (m < 1 || n < 1){
		printf ("Dimensiones incorrectas, vuelva a teclear las dimensiones\n");
		scanf ("%d %d", &m, &n);
	}
	printf ("%d %d\n", m, n);
	
	/* ---------------------------------------------------------------------- */
	
	/* RESERVA DE MEMORIA */
	x = (double *) malloc (n * sizeof (double));
	y = (double *) malloc (n * sizeof (double));
	z = (double *) malloc (n * sizeof (double));
	b = (double *) malloc (m * sizeof (double));
    Ax = (double *) malloc (m * sizeof(double));
    Atb = (double *) malloc (n * sizeof (double));
	A = (double **) malloc (m * sizeof (double *));
	At = (double **) malloc (n * sizeof (double *));
    AtA = (double **) malloc (n * sizeof (double *));
	if (x == NULL || y == NULL || z == NULL || b == NULL || Atb == NULL ||
		A == NULL || At == NULL || AtA == NULL){
		
		printf ("No hay suficiente memoria\n");
		return 1;
	}
	for (i = 0; i < m; i++){
		A[i] = (double *) malloc (n * sizeof (double));
		if (A[i] == NULL){
			printf ("No hay suficiente memoria\n");
			return 1;
		}
	}
	for (i = 0; i < n; i++){
		At[i] = (double *) malloc (m * sizeof (double));
		if (At[i] == NULL){
			printf ("No hay suficiente memoria\n");
			return 1;
		}
	}
	for (i = 0; i < n; i++){
		AtA[i] = (double *) malloc (n * sizeof (double));
		if (AtA[i] == NULL){
			printf ("No hay suficiente memoria\n");
			return 1;
		}
	}
	
	/* ---------------------------------------------------------------------- */
	
	/* LECTURA DE DATOS */
	printf ("\nMatriz A:\n");
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			scanf ("%le", &A[i][j]);
			printf ("%le ", A[i][j]);
		}
		printf ("\n");
	}
	printf ("\nMatriu A traspuesta\n");
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			At[i][j] = A[j][i];	
            printf ("%le ", At[i][j]);
		}
		printf ("\n");
	}
	printf ("\nVector b:\n");
	for (i = 0; i < m; i++){
		scanf ("%le", &b[i]);
		printf ("%le\n", b[i]);
	}
	printf ("\nTolerancia:");
	scanf ("%le", &tolerancia);
	printf ("%le\n", tolerancia);
	
	/* ---------------------------------------------------------------------- */
	
	/* CALCULO DE LA SOLUCION */
	prodMatMat (n, m, n, At, A, AtA);
	/* Multiplicamos por la izquierda la matriz A, para transformarla en
	 * cuadrada */
	prodMatVec (n, m, At, b, Atb);
	/* Por lo tanto, hemos de hacer lo mismo con el otro lado del sistema */
	if (ldlt (n, AtA, tolerancia) == 1){
		/* Hacemos la descomposicion, si la funcion ldlt retorna 1, es que no se
		 * ha podido realizar */
		printf ("No se pudo descomponer\n");
		return 0;
	}
	resTinf (n, AtA, z, Atb);	/* Resolvemos el sistema Lz = b */
	for (i = 0; i < n; i++){	/* Resolvemos el sistema Uy = z */
		y[i] = z[i] / AtA[i][i];
	}
	resTsup (n, AtA, x, y);		/* Resolvemos el sistema Ltz = x */
	printf ("Solucion\n");
	/* Imprimimos la solucion */
	for (i = 0; i < m; i++){
        printf ("x%d = ", i);
		printf ("%le\n", x[i]);
    }
    
    /* ---------------------------------------------------------------------- */
	
	/* CALCULO DE LA NORMA */
	prodMatVec (m, n, A, x, Ax);
	/* Calculamos Ax */
	for (i = 0; i < m; i++){
		norma += (Ax[i]-b[i]) * (Ax[i]-b[i]);
		/* Para cada termino sumamos a la norma (Ax-b)² */
	}
	printf ("\n\nNorma matricial: %23.16le\n", sqrt(norma));
	/* Imprimimos la norma, despues de hacer la raiz cuadrada */
	
	/* ---------------------------------------------------------------------- */
	
	/* LIBERACION DE MEMORIA */
    free (Ax);
	free (z);
	free (y);
	free (x);
	free (b);
	free (Atb);
	for (i = 0; i < n; i++){
		free (AtA[i]);
		free (At[i]);
	}
	free (At);
	free (AtA);
	for (i = 0; i < m; i++){
		free (A[i]);
	}
	free (A);
	return 0;
}

int ldlt (int n, double **A, double tol){
	
	/* DECLARACION DE VARIABLES */
	int i, j, k;
	/* Varaibles auxiliares */
	double sum;
	/* Aqui guardaremos el valor de los sumatorios */
	
	/* CALCULO DESCOMPOSICION */
	for (k = 0; k < n; k++){
		
		sum = 0;
		
		/* Calculo termino diagonal segun la formula */
		for (j = 0; j < k; j++){
			sum += A[k][j]*A[k][j]*A[j][j];
		}
		A[k][k] -= sum;
		if (fabs(A[k][k]) < tol){
			/* Si el termino diagonal es < que una tolerancia, salimos, no se
			 * pudo descomponer */
			return 1;
		}
		
		for (i = k+1; i < n; i++){
			sum = 0;
			
			/* Calculo de cada elemento de la matriz inferior L segun la 
			 * formula */
			for (j = 0; j < k; j++){
				sum += A[i][j]*A[k][j]*A[j][j];
			}
			A[i][k] -= sum;
			A[i][k] /= A[k][k];
            A[k][i] = A[i][k];
			/* Generamos Lt segun la definicion de traspuesta a partir de L */
		}
	}
	return 0;
}
