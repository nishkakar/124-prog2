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

void print_graph(float** matrix, int dimension) {
    // prints out adjacency matrix
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            printf("%.0f ", matrix[i][j]);
        }
        printf("\n");
    }
}

float** strassen(float** matrix_A, float** matrix_B, int dimension, int crossover_dimension) {
    return matrix_A;
}

int main(int argc, char *argv[]) {
    int dimension = atoi(argv[2]);
    char* inputfile = argv[3];

    FILE* fp;
    fp = fopen(inputfile, "r");
    float** matrix_A = construct_matrix(dimension, fp);
    float** matrix_B = construct_matrix(dimension, fp);

    // times the calculation for all possible crossover points
    time_t t;
    int crossover_dimension = 1;
    int optimal_dimension = -1;
    int best_time = 10E6;
    int total_time;
    while (crossover_dimension < dimension) {
        clock_t start = clock();
        float** product_matrix = strassen(matrix_A, matrix_B, dimension, crossover_dimension);
        total_time = (float) (clock() - start) / CLOCKS_PER_SEC;

        if (crossover_dimension == 1) {
            print_graph(product_matrix, dimension);
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