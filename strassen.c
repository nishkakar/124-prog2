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


float** construct_matrix(int dimension, FILE* fp) {
    float** matrix = malloc(dimension * sizeof(int*));
    char buf[256];
    for (int i = 0; i < dimension; i++) {
        matrix[i] = malloc(dimension * sizeof(int));
        for (int j = 0; j < dimension; j++) {
            fgets(buf, sizeof(buf), fp);
            matrix[i][j] = atoi(buf);
            printf("%.0f ", matrix[i][j]);
        }
        printf("\n");
    }
    return matrix;
}

// void print_graph(float** matrix, int dimension) {
//     // prints out adjacency matrix
//     for (int i = 0; i < dimension; ++i) {
//         for (int j = 0; j < dimension; ++j) {
//             printf("%f", matrix[i][j]);
//         }
//         printf("\n");
//     }
// }

int main(int argc, char *argv[]) {
    int dimension = atoi(argv[2]);
    char* inputfile = argv[3];

    FILE* fp;
    fp = fopen(inputfile, "r");
    float** matrix_A = construct_matrix(dimension, fp);
    float** matrix_B = construct_matrix(dimension, fp);

    return 0;
}