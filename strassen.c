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


matrix construct_matrix(int dimension, FILE* fp) {
    matrix m; 
    int** matrix = malloc(dimension * sizeof(int*));

    char buf[256];
    for (int i = 0; i < dimension; i++) {
        matrix[i] = malloc(dimension * sizeof(int));
        for (int j = 0; j < dimension; j++) {
            fgets(buf, sizeof(buf), fp);
            matrix[i][j] = atoi(buf);
            // printf("%d ", matrix[i][j]);
        }
        // printf("\n");
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
    for (int i = fr; i < lr; ++i) {
        for (int j = fc; j < lc; ++j) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

// int add(int** matrix_A, int** matrix_B, int dimension) {
//     int sum_matrix[dimension][dimension];
//     for (int i = 0; i < dimension; i++) {
//         for (int j = 0; j < dimension; j++) {
//             sum_matrix[i][j] = matrix_A[i][j] + matrix_B[i][j];
//         }
//     }
//     return sum_matrix;
// }

matrix strassen(matrix A, matrix B, int dimension, int crossover_dimension) {
    if (dimension < crossover_dimension) {
        // return standard_multiplication(matrix_A, matrix_B);
    }

    int** array[9];
    // array[0] = F-H; // diff(F, H, array[0])
    // array[0] = strassen(A, array[0]); // P1
    // array[1] = A+B;
    // array[1] = strassen(array[1], H); // P2
    // array[2] = C+D;
    // array[2] = strassen(array[2], E); // P3
    // array[3] = G-E;
    // array[3] = strassen(D, array[3]); // P4
    // array[4] = A+D;
    // array[5] = E+H;
    // array[4] = strassen(array[4], array[5]); // P5
    // array[5] = B-D;
    // array[6] = G+H;
    // array[5] = strassen(array[5], array[6]); // P6
    // array[6] = A-C;
    // array[7] = E+F;
    // array[6] = strassen(array[6], array[7]); // P7
    // array[7] = array[4] + array[3] - array[1] + array[5]; // AE + BG
    // array[5] = array[0] + array[1]; // AF + BH
    // array[1] = array[2] + array[3]; // CE + DG
    // array[3] = array[4] + array[0] - array[2] - array[6] // CF + DH

    return A;
}

int main(int argc, char *argv[]) {
    int dimension = atoi(argv[2]);
    char* inputfile = argv[3];

    FILE* fp;
    fp = fopen(inputfile, "r");
    matrix A = construct_matrix(dimension, fp);
    matrix B = construct_matrix(dimension, fp);

    print_graph(A.mat, dimension);

    // times the calculation for all possible crossover points
    time_t t;
    int crossover_dimension = 1;
    int optimal_dimension = -1;
    int best_time = 10E6;
    int total_time;
    while (crossover_dimension < dimension) {
        clock_t start = clock();
        matrix product_matrix = strassen(A, B, dimension, crossover_dimension);
        total_time = (float) (clock() - start) / CLOCKS_PER_SEC;

        if (crossover_dimension == 1) {
            print_graph(product_matrix.mat, dimension);
        }

        // if this was a faster calculation, update our records
        if (total_time < best_time) {
            optimal_dimension = crossover_dimension;
            best_time = total_time;
        }
        crossover_dimension++;
    }

    return 0;
}