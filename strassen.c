// ===============================
// CS124 Programming Assignment 2
// Ajay Nathan & Nishant Kakar
// ===============================

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef struct matrix {
    int fr; // the first row we want
    int lr; // the last row we don't want
    int fc; // the first column we want
    int lc; // the last column we don't want
    int** mat;
} matrix;

// frees the matrix of a given dimension
void free_matrix(int** matrix, int dimension) {
    for (int i = 0; i < dimension; ++i) {
        free(matrix[i]);
    }

    free(matrix);
}

// initializes matrix of a given dimension
int** initialize_matrix(int dimension) {
    int** mat = malloc(dimension * sizeof(int*));

    for (int i = 0; i < dimension; i++) {
        mat[i] = malloc(dimension * sizeof(int));
    }

    return mat;
}

// checks if the given integer is a power of two
int isPowerOfTwo(int n){
    if (n == 0) 
        return 0;
    while (n != 1) {
        if (n % 2 != 0)
            return 0;
        n = n/2;
    }
    return 1;
}

// true dimension is the dimension with padding; dimension is the actual dimension of the matrix we care about
// constructs matrix with the passed in file and dimensions
matrix construct_matrix(int dimension, int true_dimension, FILE* fp) {
    matrix m; 
    
    int** matrix = initialize_matrix(true_dimension);

    char buf[256];
    for (int i = 0; i < true_dimension; i++) {
        for (int j = 0; j < true_dimension; j++) {

            if (i >= dimension || j >= dimension) {
                matrix[i][j] = 0;
                continue;
            }
            fgets(buf, sizeof(buf), fp);
            matrix[i][j] = atoi(buf);
         }
    }

    m.fr = 0;
    m.lr = true_dimension;
    m.fc = 0;
    m.lc = true_dimension;
    m.mat = matrix;

    return m;
}

// largely for debugging purposes, prints the matrix
void print_matrix(matrix* M) {
    // prints out adjacency matrix
    for (int i = M->fr; i < M->lr; ++i) {
        for (int j = M->fc; j < M->lc; ++j) {
            printf("%d ", M->mat[i][j]);
        }
        printf("\n");
    }
}

// prints the diagonals of the matrix, ignoring the padding (uses dimension, not true_dimension)
void print_diagonals(matrix* M, int dimension) {
    // prints out diagonals
    for (int i = M->fr; i < M->fr + dimension; ++i) {
        printf("%d\n", M->mat[i][i]);
    }
}

// sums two matrices (A+B) and stores the result in C
void sum(matrix* A, matrix* B, matrix* C) {
    int dim = A->lr - A->fr;

    for (int Ai = A->fr, Bi = B->fr, i = 0; i < dim; ++Ai, ++Bi, ++i) {
        for (int Aj = A->fc, Bj = B->fc, j = 0; j < dim; ++Aj, ++Bj, ++j) {
            C->mat[i][j] = A->mat[Ai][Aj] + B->mat[Bi][Bj];
        }
    }

    C->fr = C->fc = 0;
    C->lr = C->lc = dim;
}

// subtracts two matrices (A-B) and stores the result in C
void diff(matrix* A, matrix* B, matrix* C) {
    int dim = A->lr - A->fr;

    for (int Ai = A->fr, Bi = B->fr, i = 0; i < dim; ++Ai, ++Bi, ++i) {
        for (int Aj = A->fc, Bj = B->fc, j = 0; j < dim; ++Aj, ++Bj, ++j) {
            C->mat[i][j] = A->mat[Ai][Aj] - B->mat[Bi][Bj];
        }
    }

    C->fr = C->fc = 0;
    C->lr = C->lc = dim;
}

// multiplies two matrices (A*B) and stores the result in C
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

// initializes a matrix struct with the passed in row & column values and matrix
void set_matrix(matrix* M, int fr, int lr, int fc, int lc, int** mat) {
    M->fr = fr;
    M->lr = lr;
    M->fc = fc;
    M->lc = lc;
    M->mat = mat;
}

// recursive implementation of strassen's algorithm
// multiples two matrices (M1*M2) and stores the result in 'result'
// if the dimension is ever lower than the crossover_dimension, then it will stop the recursion and perform standard multiplication
void strassen(matrix* M1, matrix* M2, matrix* result, int dimension, int crossover_dimension) {
    if (dimension <= crossover_dimension) {
        standard_multiplication(M1, M2, result);
        return;
    }

    // divides M1 and M2 into the sub-matrices according to the lecture notes
    matrix A, B, C, D, E, F, G, H;
    set_matrix(&A, 0, dimension/2, 0, dimension/2, M1->mat);
    set_matrix(&B, 0, dimension/2, dimension/2, dimension, M1->mat);
    set_matrix(&C, dimension/2, dimension, 0, dimension/2, M1->mat);
    set_matrix(&D, dimension/2, dimension, dimension/2, dimension, M1->mat);
    set_matrix(&E, 0, dimension/2, 0, dimension/2, M2->mat);
    set_matrix(&F, 0, dimension/2, dimension/2, dimension, M2->mat);
    set_matrix(&G, dimension/2, dimension, 0, dimension/2, M2->mat);
    set_matrix(&H, dimension/2, dimension, dimension/2, dimension, M2->mat);

    // initializes 9 temporary matrices that we'll use to compute P1-P7 and the final 4 quadrant matrices (which combine to form 'result')
    matrix temp_matrices[9];
    for (int i = 0; i < 9; ++i) {
        set_matrix(&temp_matrices[i], 0, dimension/2, 0, dimension/2, initialize_matrix(dimension/2));
    }

    // F-H
    diff(&F, &H, &temp_matrices[0]);
    // P1 = A(F-H)
    strassen(&A, &temp_matrices[0], &temp_matrices[1], dimension/2, crossover_dimension);

    // A+B
    sum(&A, &B, &temp_matrices[0]);
    // P2 = (A+B)H
    strassen(&temp_matrices[0], &H, &temp_matrices[2], dimension/2, crossover_dimension);

    // C+D 
    sum(&C, &D, &temp_matrices[0]);
    // P3 = (C+D)E
    strassen(&temp_matrices[0], &E, &temp_matrices[3], dimension/2, crossover_dimension);

    // G-E  
    diff(&G, &E, &temp_matrices[0]);
    // P4 = D(G-E)
    strassen(&D, &temp_matrices[0], &temp_matrices[4], dimension/2, crossover_dimension);

    // A+D   
    sum(&A, &D, &temp_matrices[0]);
    // E+H    
    sum(&E, &H, &temp_matrices[8]);
    // P5 = (A+D)(E+H)
    strassen(&temp_matrices[0], &temp_matrices[8], &temp_matrices[5], dimension/2, crossover_dimension);

    // B-D
    diff(&B, &D, &temp_matrices[0]);
    // G+H
    sum(&G, &H, &temp_matrices[8]);
    // P6 = (B-D)(G+H)
    strassen(&temp_matrices[0], &temp_matrices[8], &temp_matrices[6], dimension/2, crossover_dimension);

    // A-C
    diff(&A, &C, &temp_matrices[0]);
    // E+F
    sum(&E, &F, &temp_matrices[8]);
    // P7 = (A-C)(E+F)
    strassen(&temp_matrices[0], &temp_matrices[8], &temp_matrices[7], dimension/2, crossover_dimension);

    // AE + BG = P5 + P4 - P2 + P6
    sum(&temp_matrices[5], &temp_matrices[4], &temp_matrices[0]);
    diff(&temp_matrices[0], &temp_matrices[2], &temp_matrices[8]);
    sum(&temp_matrices[8], &temp_matrices[6], &temp_matrices[0]);

    // AF + BH = P1 + P2
    sum(&temp_matrices[1], &temp_matrices[2], &temp_matrices[6]);

    // CE + DG = P3 + P4
    sum(&temp_matrices[3], &temp_matrices[4], &temp_matrices[2]);

    // CF + DH = P5 + P1 - P3 - P7
    sum(&temp_matrices[5], &temp_matrices[1], &temp_matrices[4]);
    diff(&temp_matrices[4], &temp_matrices[3], &temp_matrices[5]);
    diff(&temp_matrices[5], &temp_matrices[7], &temp_matrices[4]);

    // free P1, P3, P5, P7 and the temporary storage matrix; P2, P4, P6 were overwritten with the quadrant matrices so we still need those
    free_matrix(temp_matrices[1].mat, dimension/2);
    free_matrix(temp_matrices[3].mat, dimension/2);
    free_matrix(temp_matrices[5].mat, dimension/2);
    free_matrix(temp_matrices[7].mat, dimension/2);
    free_matrix(temp_matrices[8].mat, dimension/2);

    // combine quadrant matrices into 'result'
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            if (i < dimension/2 && j < dimension/2) {
                result->mat[i][j] = temp_matrices[0].mat[i][j];  
            }
            else if (i < dimension/2 && j >= dimension/2) {
                result->mat[i][j] = temp_matrices[6].mat[i][j % (dimension/2)];
            }
            else if (i >= dimension/2 && j < dimension/2) {
                result->mat[i][j] = temp_matrices[2].mat[i % (dimension/2)][j];
            }
            else {
                result->mat[i][j] = temp_matrices[4].mat[i % (dimension/2)][j % (dimension/2)];
            }
        }
    }

    // freeing the 4 quadrant storage matrices
    free_matrix(temp_matrices[0].mat, dimension/2);
    free_matrix(temp_matrices[2].mat, dimension/2);
    free_matrix(temp_matrices[4].mat, dimension/2);
    free_matrix(temp_matrices[6].mat, dimension/2);

    return;
}

int main(int argc, char* argv[]) {
    int dimension = atoi(argv[2]);
    char* inputfile = argv[3];

    FILE* fp;
    fp = fopen(inputfile, "r");

    // checks if padding is necessary; if so, adjusts true_dimension to account for padding
    int true_dimension = dimension;
    if (!isPowerOfTwo(dimension)) {
        true_dimension = (int) pow(2.0, (double) floor(log(dimension)/log(2)) + 1);
    }
    matrix A = construct_matrix(dimension, true_dimension, fp);
    matrix B = construct_matrix(dimension, true_dimension, fp);

    int crossover_dimension = true_dimension/2;
    matrix result;
    set_matrix(&result, 0, true_dimension, 0, true_dimension, initialize_matrix(true_dimension));
    strassen(&A, &B, &result, true_dimension, crossover_dimension);
    print_diagonals(&result, dimension);
    free_matrix(result.mat, true_dimension);

    // ==============================================
    // code used to find our experimental crossover value
    // times the calculation for all possible crossover points
    // ==============================================
    // time_t t;
    // float best_time = 10E6;
    // float total_time;
    // int optimal_dimension = -1;
    // while (crossover_dimension <= true_dimension) {
    //     matrix result;
    //     set_matrix(&result, 0, true_dimension, 0, true_dimension, initialize_matrix(true_dimension));
    //     clock_t start = clock();
    //     strassen(&A, &B, &result, true_dimension, crossover_dimension);
    //     total_time = (float) (clock() - start) / CLOCKS_PER_SEC;

    //     printf("PRODUCT CROSSOVER %d %f\n", crossover_dimension, total_time);
    //     free_matrix(result.mat, true_dimension);
    //     printf("\n");

    //     // if this was a faster calculation, update our records
    //     if (total_time <= best_time) {
    //         optimal_dimension = crossover_dimension;
    //         best_time = total_time;
    //     }

    //     // only check crossover points that are powers of 2
    //     // crossover_dimension *= 2;
    //     crossover_dimension++;
    // }
    // printf("OPTIMAL DIMENSION: %d\n", optimal_dimension);

    free_matrix(A.mat, true_dimension);
    free_matrix(B.mat, true_dimension);

    return 0;
}
