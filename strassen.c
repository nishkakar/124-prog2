// ===============================
// CS124 Programming Assignment 2
// Ajay Nathan & Nishant Kakar
// ===============================

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

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

int** intialize_matrix(int dimension) {
    int** matrix = malloc(dimension * sizeof(int) + dimension * dimension * sizeof(int));
    int* pos = (int*) (matrix + dimension);
    for (int i = 0; i < dimension; i++) {
        matrix[i] = pos + i * dimension;
    }
    return matrix;
}

matrix construct_matrix(int dimension, FILE* fp) {
    matrix m; 
    int** matrix = malloc(dimension * sizeof(int*) + dimension * dimension * sizeof(int));

    int* pos = (int*) (matrix + dimension);
    for (int i = 0; i < dimension; ++i) {
        matrix[i] = pos + i * dimension;
    }

    char buf[256];
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            fgets(buf, sizeof(buf), fp);
            matrix[i][j] = atoi(buf);
         }
    }

    m.fr = 0;
    m.lr = dimension;
    m.fc = 0;
    m.lc = dimension;
    m.mat = matrix;

    return m;
}

void print_matrix(matrix M) {
    // prints out adjacency matrix
    for (int i = M.fr; i < M.lr; ++i) {
        for (int j = M.fc; j < M.lc; ++j) {
            printf("%d ", M.mat[i][j]);
        }
        printf("\n");
    }
}

matrix sum(matrix A, matrix B, matrix C) {
    printf("entered sum\n"); 

    for (int a_i = A.fr, b_i = B.fr, i = 0; a_i < A.lr; a_i++, b_i ++, i++) {
        printf("entered first for loop\n");
        for (int a_j = A.fc, b_j = B.fc, j = 0; a_j < A.lc; a_j++, b_j++, j++) {
            printf("entered second for loop\n");
            printf("i=%d, j=%d\n", i, j);
            C.mat[i][j] = A.mat[a_i][a_j] + B.mat[b_i][b_j];
        }
    }

    return C;
}

matrix diff(matrix A, matrix B, matrix C) {
    for (int a_i = A.fr, b_i = B.fr, i = 0; a_i < A.lr; a_i++, b_i++, i++) {
        for (int a_j = A.fc, b_j = B.fc, j = 0; a_j < A.lc; a_j++, b_j++, j++) {
            printf("entering %d %d\n", A.mat[a_i][a_j], B.mat[b_i][b_j]);
            C.mat[i][j] = A.mat[a_i][a_j] - B.mat[b_i][b_j];
            printf("exiting\n");
        }
    }

    return C;
}

// AB = C
void standard_multiplication(matrix A, matrix B, matrix* C) {
	int dimension = A.lr - A.fr, sum = 0;

	C->fr = 0;
	C->lr = dimension;
	C->fc = 0;
	C->lc = dimension;

	for (int i = 0, Ai = A.fr; i < dimension; ++i, ++Ai) {
		for (int j = 0, Bj = B.fc; j < dimension; ++j, ++Bj) {
			for (int Ak = A.fc, Bk = B.fr; Ak < A.lc; ++Ak, ++Bk) {
				sum += A.mat[Ai][Ak] * B.mat[Bk][Bj];
			}

			C->mat[i][j] = sum;
			sum = 0;
		}
	}
}

matrix strassen(matrix M1, matrix M2, int dimension, int crossover_dimension) {
    if (dimension <= crossover_dimension) {
        printf("STANDARD\n");
    	matrix temp_matrix = {.fr = 0, .lr = dimension, .fc = 0, .lc = dimension, intialize_matrix(dimension)};
        standard_multiplication(M1, M2, &temp_matrix);
        return temp_matrix;
    }

    matrix A = {.fr = 0, .lr = dimension/2, .fc = 0, .lc = dimension/2, .mat = M1.mat};
    matrix B = {.fr = 0, .lr = dimension/2, .fc = dimension/2, .lc = dimension, .mat = M1.mat};
    matrix C = {.fr = dimension/2, .lr = dimension, .fc = 0, .lc = dimension/2, .mat = M1.mat};
    matrix D = {.fr = dimension/2, .lr = dimension, .fc = dimension/2, .lc = dimension, .mat = M1.mat};
    matrix E = {.fr = 0, .lr = dimension/2, .fc = 0, .lc = dimension/2, .mat = M2.mat};
    matrix F = {.fr = 0, .lr = dimension/2, .fc = dimension/2, .lc = dimension, .mat = M2.mat};
    matrix G = {.fr = dimension/2, .lr = dimension, .fc = 0, .lc = dimension/2, .mat = M2.mat};
    matrix H = {.fr = dimension/2, .lr = dimension, .fc = dimension/2, .lc = dimension, .mat = M2.mat};

    matrix temp_matrices[9];
    for (int i = 0; i < 9; i++) {
        temp_matrices[i].fr = 0;
        temp_matrices[i].lr = dimension/2;
        temp_matrices[i].fc = 0;
        temp_matrices[i].lc = dimension/2; 
        temp_matrices[i].mat = intialize_matrix(dimension/2);
    }

    printf("INITIALIZED\n");

    // array[0] = F-H; // diff(F, H, array[0])
    printf("===================%p\n", temp_matrices[1].mat[0]);
    diff(F, H, temp_matrices[0]);
    printf("===================%p\n", temp_matrices[1].mat[0]);

    // array[0] = strassen(A, array[0]); // P1
    temp_matrices[0] = strassen(A, temp_matrices[0], dimension/2, crossover_dimension);
    printf("Done with first strassen\n");

    // array[1] = A+B;
    sum(A, B, temp_matrices[1]);
    printf("Done with sum\n");


    // array[1] = strassen(array[1], H); // P2
    temp_matrices[1] = strassen(temp_matrices[1], H, dimension/2, crossover_dimension);
    printf("Done with second strassen\n");


    // array[2] = C+D;
    sum(C, D, temp_matrices[2]);

    // array[2] = strassen(array[2], E); // P3
    temp_matrices[2] = strassen(temp_matrices[2], E, dimension/2, crossover_dimension);

    // array[3] = G-E;
    diff(G, E, temp_matrices[3]);

    // array[3] = strassen(D, array[3]); // P4
    temp_matrices[3] = strassen(D, temp_matrices[3], dimension/2, crossover_dimension);

    // array[4] = A+D;
    sum(A, D, temp_matrices[4]);

    // array[5] = E+H;
    sum(E, H, temp_matrices[5]);

    // array[4] = strassen(array[4], array[5]); // P5
    temp_matrices[4] = strassen(temp_matrices[4], temp_matrices[5], dimension/2, crossover_dimension);

    // array[5] = B-D;
    diff(B, D, temp_matrices[5]);

    // array[6] = G+H;
    sum(G, H, temp_matrices[6]);

    // array[5] = strassen(array[5], array[6]); // P6
    temp_matrices[5] = strassen(temp_matrices[5], temp_matrices[6], dimension/2, crossover_dimension);

    // array[6] = A-C;
    diff(A, C, temp_matrices[6]);

    // array[7] = E+F;
    sum(E, F, temp_matrices[7]);

    // array[6] = strassen(array[6], array[7]); // P7
    temp_matrices[6] = strassen(temp_matrices[6], temp_matrices[7], dimension/2, crossover_dimension);

    // array[7] = array[4] + array[3] - array[1] + array[5]; // AE + BG
    sum(diff(sum(temp_matrices[4], temp_matrices[3], temp_matrices[7]), temp_matrices[1], temp_matrices[8]), temp_matrices[5], temp_matrices[7]);

    // array[5] = array[0] + array[1]; // AF + BH
    sum(temp_matrices[0], temp_matrices[1], temp_matrices[5]);

    // array[1] = array[2] + array[3]; // CE + DG
    sum(temp_matrices[2], temp_matrices[3], temp_matrices[1]);

    // array[3] = array[4] + array[0] - array[2] - array[6] // CF + DH
    diff(diff(sum(temp_matrices[4], temp_matrices[0], temp_matrices[3]), temp_matrices[2], temp_matrices[8]), temp_matrices[6], temp_matrices[3]);

    printf("COMBINE\n");

    // combine matrices
    matrix M = {.fr = 0, .lr = dimension, .fc = 0, .lc = dimension, intialize_matrix(dimension)};
    for (int i = 0; i < dimension; ++i) {
    	for (int j = 0; j < dimension; ++j) {
    		if (i < dimension/2 && j < dimension/2) {
    			M.mat[i][j] = temp_matrices[7].mat[i][j];  			
    		}
    		else if (i < dimension/2 && j >= dimension/2) {
    			M.mat[i][j] = temp_matrices[5].mat[i][j % (dimension/2)];
    		}
    		else if (i >= dimension/2 && j < dimension/2) {
    			M.mat[i][j] = temp_matrices[1].mat[i % (dimension/2)][j];
    		}
    		else {
    			M.mat[i][j] = temp_matrices[3].mat[i % (dimension/2)][j % (dimension/2)];
    		}
    	}
    }

    printf("DONE COMBINING\n");

    return M;
}

int main(int argc, char *argv[]) {
    int dimension = atoi(argv[2]);
    char* inputfile = argv[3];

    FILE* fp;
    fp = fopen(inputfile, "r");
    matrix A = construct_matrix(dimension, fp);
    matrix B = construct_matrix(dimension, fp);

    // times the calculation for all possible crossover points
    time_t t;
    int crossover_dimension = 2;
    int optimal_dimension = -1;
    int best_time = 10E6;
    int total_time;
    while (crossover_dimension <= dimension) {
        clock_t start = clock();
        printf("ENTERING STRASSEN\n");
        matrix product_matrix = strassen(A, B, dimension, crossover_dimension);
        printf("LEAVING STRASSEN\n");
        total_time = (float) (clock() - start) / CLOCKS_PER_SEC;
	
		printf("PRODUCT CROSSOVER %d\n", crossover_dimension);
		print_matrix(product_matrix);
		printf("\n");

        // if this was a faster calculation, update our records
        if (total_time < best_time) {
            optimal_dimension = crossover_dimension;
            best_time = total_time;
        }
        crossover_dimension++;
    }

    matrix C = {.fr = 0, .lr = dimension, .fc = 0, .lc = dimension, intialize_matrix(dimension)};
    
    standard_multiplication(A, B, &C);
    printf("STANDARD MULT\n");
    print_matrix(C);
    printf("\n");

    return 0;
}
