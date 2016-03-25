// ===============================
// CS124 Programming Assignment 2
// Ajay Nathan & Nishant Kakar
// ===============================

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// #include <string.h>

// for seeking
// #include <sys/types.h>
// #include <unistd.h>

typedef struct matrix {
    int fr; // the first row we want
    int lr; // the last row we don't want
    int fc; // the first column we want
    int lc; // the last column we don't want
    int** mat;
} matrix; 

int** initialize_matrix(int dimension) {
    int** mat = malloc(dimension * sizeof(int*));

    for (int i = 0; i < dimension; i++) {
        mat[i] = malloc(dimension * sizeof(int));
    }

    return mat;
}

matrix construct_matrix(int dimension, FILE* fp) {
    matrix m; 
    int** matrix = initialize_matrix(dimension); // malloc(dimension * sizeof(int*) + dimension * dimension * sizeof(int));

    // int* pos = (int*) (matrix + dimension);
    // for (int i = 0; i < dimension; ++i) {
    //     matrix[i] = pos + i * dimension;
    // }

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

matrix construct_matrix2(int dimension, FILE* fp) {
    matrix m;
    dimension = 2 * dimension;
    int** matrix = malloc(dimension * sizeof(int*) + dimension * dimension * sizeof(int));

    int* pos = (int*) (matrix + dimension);
    for (int i = 0; i < dimension; ++i) {
        matrix[i] = pos + i * dimension;
    }

    char buf[256];
    for (int i = dimension/2; i < dimension; i++) {
        for (int j = dimension; j < dimension; j++) {
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

// int** initialize_matrix(int dimension) {
//     int** mat = malloc(dimension * sizeof(int*) + dimension * dimension * sizeof(int));

//     int* pos = (int*) (mat + dimension);
//     for (int i = 0; i < dimension * dimension; i += dimension) {
//         mat[i] = pos + i;
//     }

//     return mat;
// }

void print_matrix(matrix* M) {
    // prints out adjacency matrix
    for (int i = M->fr; i < M->lr; ++i) {
        for (int j = M->fc; j < M->lc; ++j) {
            printf("%d ", M->mat[i][j]);
        }
        printf("\n");
    }
}

void print_diagonals(matrix* M) {
    // prints out diagonals
    for (int i = M->fr; i < M->lr; ++i) {
        printf("%d\n", M->mat[i][i]);
    }
}

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

void strassen(matrix* M1, matrix* M2, matrix* result, int dimension, int crossover_dimension) {
    if (dimension <= crossover_dimension) {
        standard_multiplication(M1, M2, result);
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
    // array[1] = strassen(&A, array[0]); // P1
    strassen(&A, &temp_matrices[0], &temp_matrices[1], dimension/2, crossover_dimension);
    // array[0] = A+B;
    sum(&A, &B, &temp_matrices[0]);
    // array[2] = strassen(array[0], H); // P2
    strassen(&temp_matrices[0], &H, &temp_matrices[2], dimension/2, crossover_dimension);
    // array[0] = C+D;   
    sum(&C, &D, &temp_matrices[0]);
    // array[3] = strassen(array[0], E); // P3
    strassen(&temp_matrices[0], &E, &temp_matrices[3], dimension/2, crossover_dimension);
    // array[0] = G-E;    
    diff(&G, &E, &temp_matrices[0]);
    // array[4] = strassen(&D, array[0]); // P4
    strassen(&D, &temp_matrices[0], &temp_matrices[4], dimension/2, crossover_dimension);
    // array[0] = A+D;   
    sum(&A, &D, &temp_matrices[0]);
    // array[8] = E+H;    
    sum(&E, &H, &temp_matrices[8]);
    // array[5] = strassen(array[0], array[8]); // P5
    strassen(&temp_matrices[0], &temp_matrices[8], &temp_matrices[5], dimension/2, crossover_dimension);
    // array[0] = B-D;
    diff(&B, &D, &temp_matrices[0]);
    // array[8] = G+H;
    sum(&G, &H, &temp_matrices[8]);
    // array[6] = strassen(array[0], array[8]); // P6
    strassen(&temp_matrices[0], &temp_matrices[8], &temp_matrices[6], dimension/2, crossover_dimension);
    // array[0] = A-C;
    diff(&A, &C, &temp_matrices[0]);
    // array[8] = E+F;
    sum(&E, &F, &temp_matrices[8]);
    // array[7] = strassen(array[0], array[8]); // P7
    strassen(&temp_matrices[0], &temp_matrices[8], &temp_matrices[7], dimension/2, crossover_dimension);
    // array[0] = array[5] + array[4] - array[2] + array[6]; // AE + BG
    sum(&temp_matrices[5], &temp_matrices[4], &temp_matrices[0]);
    diff(&temp_matrices[0], &temp_matrices[2], &temp_matrices[8]);
    sum(&temp_matrices[8], &temp_matrices[6], &temp_matrices[0]);
    // array[6] = array[1] + array[2]; // AF + BH
    sum(&temp_matrices[1], &temp_matrices[2], &temp_matrices[6]);
    // array[2] = array[3] + array[4]; // CE + DG
    sum(&temp_matrices[3], &temp_matrices[4], &temp_matrices[2]);
    // array[4] = array[5] + array[1] - array[3] - array[7] // CF + DH
    sum(&temp_matrices[5], &temp_matrices[1], &temp_matrices[4]);
    diff(&temp_matrices[4], &temp_matrices[3], &temp_matrices[5]);
    diff(&temp_matrices[5], &temp_matrices[7], &temp_matrices[4]);

    // combine
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

    return;
}

int main(int argc, char* argv[]) {
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
    float best_time = 10E6;
    float total_time;
    while (crossover_dimension <= dimension) {
        matrix result;
        int** temp = initialize_matrix(dimension);
        set_matrix(&result, 0, dimension, 0, dimension, temp);
        clock_t start = clock();
        strassen(&A, &B, &result, dimension, crossover_dimension);
        total_time = (float) (clock() - start) / CLOCKS_PER_SEC;
        printf("PRODUCT CROSSOVER %d %f\n", crossover_dimension, total_time);
        print_diagonals(&result);
        free(temp);
        printf("\n");

        // if this was a faster calculation, update our records
        if (total_time <= best_time) {
            optimal_dimension = crossover_dimension;
            best_time = total_time;
        }
        crossover_dimension++;
    }

    printf("OPTIMAL DIMENSION: %d\n", optimal_dimension);

    return 0;
}

// every other recurisive strassen call fails, only works for P1 - P7; dig in to the recursive calls for the even ones and see where it fucks up
// standard multiplication - fr,fc element in A was being overwritten (to -640 from -8)