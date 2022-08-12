#include "s21_matrix.h"

// main functions
matrix_t s21_create_matrix(int rows, int columns) {
    matrix_t new_matrix;
    if (rows <= 0 || columns <= 0) {
        new_matrix.columns = 0;
        new_matrix.rows = 0;
        new_matrix.matrix_type = INCORRECT_MATRIX;
    } else {
        new_matrix.columns = columns;
        new_matrix.rows = rows;
        new_matrix.matrix_type = ZERO_MATRIX;
        new_matrix.matrix = (double**)malloc(sizeof(double*) * rows);
        for (int i = 0; i < rows; i++) {
            new_matrix.matrix[i] = (double*)malloc(sizeof(double) * columns);
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                new_matrix.matrix[i][j] = 0;
            }
        }
    }
    return new_matrix;
}

void s21_remove_matrix(matrix_t *A) {
    if (A->matrix_type != INCORRECT_MATRIX && A->rows != 0) {
        for (int i = 0; i < A->rows; i++)
            free(A->matrix[i]);
        free(A->matrix);
    }
    A->rows = 0;
    A->columns = 0;
    A->matrix_type = INCORRECT_MATRIX;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    int result = SUCCESS;
    if (A->matrix_type != INCORRECT_MATRIX &&
        B->matrix_type != INCORRECT_MATRIX &&
        A->rows == B->rows &&
        A->columns == B->columns) {
        for (int i = 0; i < A->rows && result; i++) {
            for (int j = 0; j < A->columns && result; j++) {
                if (s21_compare_double(A->matrix[i][j], B->matrix[i][j]))
                    result = FAILURE;
            }
        }
    } else {
        result = FAILURE;
    }
    return result;
}

matrix_t s21_sum_matrix(matrix_t *A, matrix_t *B) {
    matrix_t result;
    result.matrix_type = INCORRECT_MATRIX;
    if (A->matrix_type != INCORRECT_MATRIX &&
        B->matrix_type != INCORRECT_MATRIX &&
        A->rows == B->rows &&
        A->columns == B->columns) {
        result = s21_create_matrix(A->rows, A->columns);
        for (int i = 0; (i < A->rows); i++) {
            for (int j = 0; (j < A->columns); j++) {
                result.matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
            }
        }
    }
    s21_check_type_matrix(&result);
    return result;
}

matrix_t s21_sub_matrix(matrix_t *A, matrix_t *B) {
    matrix_t result;
    result.matrix_type = INCORRECT_MATRIX;
    if (A->matrix_type != INCORRECT_MATRIX &&
        B->matrix_type != INCORRECT_MATRIX &&
        A->rows == B->rows &&
        A->columns == B->columns) {
        result = s21_create_matrix(A->rows, A->columns);
        for (int i = 0; (i < A->rows); i++) {
            for (int j = 0; (j < A->columns); j++) {
                result.matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
            }
        }
    }
    s21_check_type_matrix(&result);
    return result;
}

matrix_t s21_mult_number(matrix_t *A, double number) {
    matrix_t result;
    result.matrix_type = INCORRECT_MATRIX;
    if (A->matrix_type != INCORRECT_MATRIX) {
        result = s21_create_matrix(A->rows, A->columns);
        for (int i = 0; (i < A->rows); i++) {
            for (int j = 0; (j < A->columns); j++) {
                result.matrix[i][j] = A->matrix[i][j] * number;
            }
        }
    }
    s21_check_type_matrix(&result);
    return result;
}

matrix_t s21_mult_matrix(matrix_t *A, matrix_t *B) {
    matrix_t result;
    result.matrix_type = INCORRECT_MATRIX;
    if (A->matrix_type != INCORRECT_MATRIX &&
        B->matrix_type != INCORRECT_MATRIX &&
        A->rows == B->columns &&
        A->columns == B->rows) {
        result = s21_create_matrix(A->rows, B->columns);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < B->columns; j++) {
                for (int k = 0; k < A->columns; k++) {
                    result.matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
                }
            }
        }
    }
    s21_check_type_matrix(&result);
    return result;
}

matrix_t s21_transpose(matrix_t *A) {
    matrix_t result;
    result.matrix_type = INCORRECT_MATRIX;
    if (A->matrix_type != INCORRECT_MATRIX) {
        result = s21_create_matrix(A->columns, A->rows);
        for (int i = 0; i < A->columns; i++) {
            for (int j = 0; j < A->rows; j++) {
                result.matrix[i][j] = A->matrix[j][i];
            }
        }
    }
    s21_check_type_matrix(&result);
    return result;
}

matrix_t s21_calc_complements(matrix_t *A) {
    matrix_t result;
    result.matrix_type = INCORRECT_MATRIX;
    if (A->matrix_type != INCORRECT_MATRIX &&
        A->rows == A->columns && A->columns != 0) {
        result = s21_create_matrix(A->rows, A->columns);
        if (A->rows == 1) {
            result.matrix[0][0] = 1;
        } else {
            for (int i = 0; i < A->rows; i++) {
                for (int j = 0; j < A->columns; j++) {
                    matrix_t temp = s21_get_minor_matrix(A, i, j);
                    result.matrix[i][j] = s21_determinant(&temp) * pow(-1, (i + 1) + (j + 1));
                    s21_remove_matrix(&temp);
                }
            }
        }
    }
    s21_check_type_matrix(&result);
    return result;
}

double s21_determinant(matrix_t *A) {
    double result = 0 / 0.0;
    if (A->matrix_type != INCORRECT_MATRIX &&
        A->rows == A->columns) {
        if (A->rows == 1) {
            result = A->matrix[0][0];
        } else if (A->rows == 2) {
            result = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
        } else {
            double det = 0.0;
            for (int i = 0; i < A->columns; i++) {
                matrix_t temp = s21_get_minor_matrix(A, 0, i);
                double k = pow(-1, i + 2);
                det += k * A->matrix[0][i] * s21_determinant(&temp);
                s21_remove_matrix(&temp);
            }
            result = det;
        }
    }
    return result;
}

matrix_t s21_inverse_matrix(matrix_t *A) {
    matrix_t result;
    result.matrix_type = INCORRECT_MATRIX;
    if (A->matrix_type != INCORRECT_MATRIX &&
        A->rows == A->columns) {
        double dif = s21_determinant(A);
        if (s21_compare_double(0.0, dif)) {
            matrix_t trans = s21_transpose(A);
            matrix_t temp = s21_calc_complements(&trans);
            result = s21_mult_number(&temp, (1 / dif));
            s21_remove_matrix(&trans);
            s21_remove_matrix(&temp);
        } else {
            s21_remove_matrix(&result);
        }
    }
    s21_check_type_matrix(&result);
    return result;
}

// help functions
int s21_compare_double(double first, double second) {
    int result = 0;
    double dif = fabs(first - second);
    if (dif > E) result = 1;
    return result;
}

void s21_check_type_matrix(matrix_t *A) {
    if (A->matrix_type != INCORRECT_MATRIX) {
        matrix_t temp_mat = s21_create_matrix(A->rows, A->columns);
        if (s21_eq_matrix(&temp_mat, A)) {
            A->matrix_type = ZERO_MATRIX;
        } else if (A->rows == A->columns) {
            for (int i = 0; i < A->rows; i++) {
                temp_mat.matrix[i][i] = 1;
            }
            if (s21_eq_matrix(&temp_mat, A)) {
                A->matrix_type = IDENITY_MATRIX;
            }
        } else {
            A->matrix_type = CORRECT_MATRIX;
        }
        s21_remove_matrix(&temp_mat);
    }
}

matrix_t s21_get_minor_matrix(matrix_t *A, int index, int jndex) {
    matrix_t result;
    result.matrix_type = INCORRECT_MATRIX;
    if (A->matrix_type != INCORRECT_MATRIX &&
        A->rows == A->columns) {
        result = s21_create_matrix(A->rows - 1, A->columns - 1);
        for (int i = 0, mi = 0; i < A->rows; i++) {
            for (int j = 0, mj = 0; j < A->columns; j++) {
                if (i != index && j != jndex) result.matrix[mi][mj++] = A->matrix[i][j];
            }
            if (index != i) mi++;
        }
    }
    s21_check_type_matrix(&result);
    return result;
}
