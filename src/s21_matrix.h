#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SUCCESS 1
#define FAILURE 0
#define E 1e-8

typedef enum {
    CORRECT_MATRIX = 0,
    INCORRECT_MATRIX,
    IDENITY_MATRIX,
    ZERO_MATRIX
} matrix_type_t;

typedef struct matrix_struct {
    double** matrix;
    int rows;
    int columns;
    matrix_type_t matrix_type;
} matrix_t;

// main functions
matrix_t s21_create_matrix(int rows, int columns);      // tested
void s21_remove_matrix(matrix_t *A);                    // tested
int s21_eq_matrix(matrix_t *A, matrix_t *B);            // tested
matrix_t s21_sum_matrix(matrix_t *A, matrix_t *B);      // tested
matrix_t s21_sub_matrix(matrix_t *A, matrix_t *B);      // tested
matrix_t s21_mult_number(matrix_t *A, double number);    // tested
matrix_t s21_mult_matrix(matrix_t *A, matrix_t *B);      // tested
matrix_t s21_transpose(matrix_t *A);                    // tested
double s21_determinant(matrix_t *A);                    // tested
matrix_t s21_calc_complements(matrix_t *A);             // tested
matrix_t s21_inverse_matrix(matrix_t *A);

// help functions
int s21_compare_double(double first, double second);
void s21_check_type_matrix(matrix_t *A);
matrix_t s21_get_minor_matrix(matrix_t *A, int index, int jndex);

#endif  // SRC_S21_MATRIX_H_
