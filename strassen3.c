// ===============================
// CS124 Programming Assignment 2
// Ajay Nathan & Nishant Kakar
// ===============================

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// for seeking
#include <sys/types.h>
#include <unistd.h>

typedef struct matrix {
    int fr; // the first row we want
    int lr; // the last row we don't want
    int fc; // the first column we want
    int lc; // the last column we don't want
    int** mat;
} matrix; 

matrix* construct_matrix(int dimension, FILE* fp) {
    matrix* m; 
    int** mat = malloc(dimension * sizeof(int*) + dimension * dimension * sizeof(int));

    int* pos = (int*) (mat + dimension * sizeof(int));
    for (int i = 0; i < dimension; ++i) {
        mat[i] = pos + i * dimension;
    }

    char buf[256];
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            fgets(buf, sizeof(buf), fp);
            mat[i][j] = atoi(buf);
         }
    }

    m->fr = 0;
    m->lr = dimension;
    m->fc = 0;
    m->lc = dimension;
    m->mat = mat;

    return m;
}

int** initialize_matrix(int dimension) {
    int** mat = malloc(dimension * sizeof(int) + dimension * dimension * sizeof(int));

    int* pos = (int*) (mat + dimension * sizeof(int));
    for (int i = 0; i < dimension; i++) {
        mat[i] = pos + i * dimension;
    }

    return mat;
}

void print_matrix(matrix* M) {
    // prints out adjacency matrix
    for (int i = M->fr; i < M->lr; ++i) {
        for (int j = M->fc; j < M->lc; ++j) {
            printf("%d ", M->mat[i][j]);
        }
        printf("\n");
    }
}

void sum(matrix* A, matrix* B, matrix* C) {
	printf("%p %p %p\n", A->mat[0], B->mat[0], C->mat[0]);

	int dim = A->lr - A->fr;
	C->fr = C->fc = 0;
	C->lr = C->lc = dim;

	printf("%p %p %p\n", A->mat[0], B->mat[0], C->mat[0]);

	for (int Ai = A->fr, Bi = B->fr, i = 0; i < dim; ++Ai, ++Bi, ++i) {
		for (int Aj = A->fc, Bj = B->fc, j = 0; j < dim; ++Aj, ++Bj, ++j) {
			printf("%p %p %p\n", A->mat[0], B->mat[0], C->mat[0]);
			C->mat[i][j] = A->mat[Ai][Aj] + B->mat[Bi][Bj];
		}
	}

	printf("%p %p %p\n", A->mat[0], B->mat[0], C->mat[0]);
}

void diff(matrix* A, matrix* B, matrix* C) {
	int dim = A->lr - A->fr;
	C->fr = C->fc = 0;
	C->lr = C->lc = dim;

	for (int Ai = A->fr, Bi = B->fr, i = 0; i < dim; ++Ai, ++Bi, ++i) {
		for (int Aj = A->fc, Bj = B->fc, j = 0; j < dim; ++Aj, ++Bj, ++j) {
			C->mat[i][j] = A->mat[Ai][Aj] - B->mat[Bi][Bj];
		}
	}
}

void standard_multiplication(matrix* A, matrix* B, matrix* C) {
	int dim = A->lr - A->fr, sum = 0;

	C->fr = C->fc = 0;
	C->lr = C->lc = dim;

	for (int i = 0, Ai = A->fr; i < dim; ++i, ++Ai) {
		for (int j = 0, Bj = B->fc; j < dim; ++j, ++Bj) {
			for (int Ak = A->fc, Bk = B->fr; Ak < A->lc; ++Ak, ++Bk) {
				sum += A->mat[Ai][Ak] * B->mat[Bk][Bj];
			}

			C->mat[i][j] = sum;
			sum = 0;
		}
	}
}

void set_matrix(matrix* M, int fr, int lr, int fc, int lc, int** mat) {
	M->fr = fr;
	M->lr = lr;
	M->fc = fc;
	M->lc = lc;
	M->mat = mat;
}

void strassen(matrix* M1, matrix* M2, matrix* result, int dimension) {
	// base case
	if (dimension == 1) {
		set_matrix(result, 0, 1, 0, 1, result->mat);
		result->mat[0][0] = M1->mat[M1->fr][M1->fc - 1] * M2->mat[M2->fr][M2->fc - 1]; 
		printf("base case done\n");
		return;
	}

	matrix A, B, C, D, E, F, G, H;
    set_matrix(&A, 0, dimension/2, 0, dimension/2, M1->mat);
    set_matrix(&B, 0, dimension/2, dimension/2, dimension, M1->mat);
    set_matrix(&C, dimension/2, dimension, 0, dimension/2, M1->mat);
    set_matrix(&D, dimension/2, dimension, dimension/2, dimension, M1->mat);
    set_matrix(&E, 0, dimension/2, 0, dimension/2, M2->mat);
    set_matrix(&F, 0, dimension/2, dimension/2, dimension, M2->mat);
    set_matrix(&G, dimension/2, dimension, 0, dimension/2, M2->mat);
    set_matrix(&H, dimension/2, dimension, dimension/2, dimension, M2->mat);

    matrix temp_matrices[9];
    for (int i = 0; i < 9; ++i) {
    	set_matrix(&temp_matrices[i], 0, dimension/2, 0, dimension/2, initialize_matrix(dimension/2));
    }

    // array[0] = F-H; // diff(&F, &H, array[0])
    diff(&F, &H, &temp_matrices[0]);
    if (dimension == 8)
    	printf("%p\n", &temp_matrices[1]);
    // printf("first diff %d\n", dimension);
    // array[0] = strassen(&A, array[0]); // P1
    strassen(&A, &temp_matrices[0], &temp_matrices[0], dimension/2);
    if (dimension == 8) {
    	printf("%p\n", &temp_matrices[1]);
    	printf("exiting strassen 1\n");
	}
    // array[1] = A+B;
    sum(&A, &B, &temp_matrices[1]);
    if (dimension == 8) {
    	printf("entering strassen 2\n");
	}
    // array[1] = strassen(array[1], H); // P2
    strassen(&temp_matrices[1], &H, &temp_matrices[1], dimension/2);
    // array[2] = C+D;
    sum(&C, &D, &temp_matrices[2]);
    // array[2] = strassen(array[2], E); // P3
    strassen(&temp_matrices[2], &E, &temp_matrices[2], dimension/2);
    // array[3] = G-E;
    diff(&G, &E, &temp_matrices[3]);
    // array[3] = strassen(&D, array[3]); // P4
    strassen(&D, &temp_matrices[3], &temp_matrices[3], dimension/2);
    // array[4] = A+D;
    sum(&A, &D, &temp_matrices[4]);
    // array[5] = E+H;
    sum(&E, &H, &temp_matrices[5]);
    // array[4] = strassen(array[4], array[5]); // P5
    strassen(&temp_matrices[4], &temp_matrices[5], &temp_matrices[4], dimension/2);
    // array[5] = B-D;
    diff(&B, &D, &temp_matrices[5]);
    // array[6] = G+H;
    sum(&G, &H, &temp_matrices[6]);
    // array[5] = strassen(array[5], array[6]); // P6
    strassen(&temp_matrices[5], &temp_matrices[6], &temp_matrices[5], dimension/2);
    // array[6] = A-C;
    diff(&A, &C, &temp_matrices[6]);
    // array[7] = E+F;
    sum(&E, &F, &temp_matrices[7]);
    // array[6] = strassen(array[6], array[7]); // P7
    strassen(&temp_matrices[6], &temp_matrices[7], &temp_matrices[6], dimension/2);
    // array[7] = array[4] + array[3] - array[1] + array[5]; // AE + BG
    sum(&temp_matrices[4], &temp_matrices[3], &temp_matrices[7]);
    diff(&temp_matrices[7], &temp_matrices[1], &temp_matrices[7]);
    sum(&temp_matrices[7], &temp_matrices[5], &temp_matrices[7]);
    // array[5] = array[0] + array[1]; // AF + BH
    sum(&temp_matrices[0], &temp_matrices[1], &temp_matrices[5]);
    // array[1] = array[2] + array[3]; // CE + DG
    sum(&temp_matrices[2], &temp_matrices[3], &temp_matrices[1]);
    // array[3] = array[4] + array[0] - array[2] - array[6] // CF + DH
    sum(&temp_matrices[4], &temp_matrices[0], &temp_matrices[3]);
    diff(&temp_matrices[3], &temp_matrices[2], &temp_matrices[3]);
    diff(&temp_matrices[3], &temp_matrices[6], &temp_matrices[3]);

    printf("combining %d\n", dimension);

    // combine
    for (int i = 0; i < dimension; ++i) {
    	for (int j = 0; j < dimension; ++j) {
    		if (i < dimension/2 && j < dimension/2) {
    			result->mat[i][j] = temp_matrices[7].mat[i][j];  			
    		}
    		else if (i < dimension/2 && j >= dimension/2) {
    			result->mat[i][j] = temp_matrices[5].mat[i][j % (dimension/2)];
    		}
    		else if (i >= dimension/2 && j < dimension/2) {
    			result->mat[i][j] = temp_matrices[1].mat[i % (dimension/2)][j];
    		}
    		else {
    			result->mat[i][j] = temp_matrices[3].mat[i % (dimension/2)][j % (dimension/2)];
    		}
    	}
    }

    printf("done combining %d\n", dimension);

	return;
}

int main(int argc, char* argv[]) {
	int dimension = atoi(argv[2]);
    char* inputfile = argv[3];

    FILE* fp;
    fp = fopen(inputfile, "r");
    matrix* A = construct_matrix(dimension, fp);
    matrix* B = construct_matrix(dimension, fp);
    matrix result;
    set_matrix(&result, 0, dimension, 0, dimension, initialize_matrix(dimension));

    strassen(A, B, &result, dimension);

    print_matrix(A);
    print_matrix(B);
    print_matrix(&result);

  	return 0;
}
