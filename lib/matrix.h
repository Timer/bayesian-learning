#include <stdbool.h>
#include "list.h"

#ifndef MATRIX_H
#define MATRIX_H

typedef struct {
  int rows, cols;
  void **data;
} Matrix;

Matrix *matrix_zeros(int rows, int cols);
Matrix *matrix_create(int rows, int cols, int value);
Matrix *matrix_double_zeros(int rows, int cols);
Matrix *matrix_double_create(int rows, int cols, double value);
Matrix *matrix_range(int from, int to);
Matrix *matrix_from_file(char *file, bool transposed);
void matrix_to_file(Matrix *m, char *file);
void *matrix_element(Matrix *matrix, int row, int col);
void *matrix_element_by_index(Matrix *matrix, int index);
void matrix_set(Matrix *matrix, int value);
int matrix_prod(Matrix *matrix);
int matrix_sum(Matrix *matrix);
double matrix_double_sum(Matrix *matrix);
void matrix_double_set(Matrix *matrix, double value);
Matrix *matrix_sub_indices(Matrix *o_matrix, int row_start, int row_end, int col_start, int col_end);
Matrix *matrix_sub_lists(Matrix *o_matrix, List *rows, List *cols);
Matrix *matrix_sub_list_index(Matrix *o_matrix, List *rows, int col_start, int col_end);
Matrix *matrix_sub_index_list(Matrix *o_matrix, int row_start, int row_end, List *cols);
Matrix *matrix_sub_col(Matrix *matrix, int col);
Matrix *matrix_sub_row(Matrix *matrix, int row);
Matrix *matrix_sub_concat_rows(Matrix *matrix, Matrix *rows, bool parallel);
List *matrix_find_by_value(Matrix *matrix, int value);
List *matrix_double_find_by_value(Matrix *matrix, double value);
List *matrix_to_list(Matrix *matrix);
List *matrix_double_to_list(Matrix *matrix);
Matrix *matrix_from_array(int rows, int cols, int arr[]);
Matrix *matrix_from_list(List *l);
void matrix_delete(Matrix *matrix);
void matrix_scrap(Matrix *matrix);
void *matrix_element_n_dim(Matrix *m, Matrix *ind, Matrix *dims);
void matrix_display(Matrix *matrix);
void matrix_double_display(Matrix *matrix);
Matrix *matrix_add(Matrix *a, Matrix *b);
void matrix_add_in(Matrix *a, Matrix *b);
Matrix *matrix_double_add(Matrix *a, Matrix *b);
Matrix *matrix_add_int_double(Matrix *a, Matrix *b);
Matrix *matrix_sum_n_cols(Matrix *matrix, int cols);
Matrix *matrix_double_sum_n_cols(Matrix *matrix, int cols);
Matrix *matrix_double_subtract(Matrix *m1, Matrix *m2);
Matrix *matrix_lgamma(Matrix *m);
Matrix *matrix_create_sz(Matrix *matrix);
void matrix_double_mk_stochastic(Matrix *m, Matrix *ns);

#endif
