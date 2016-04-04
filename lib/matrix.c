#include "matrix.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// --- PRIVATE START
int _matrix_index_for(Matrix *m, int row, int col) { return col * m->rows + row; }
// --- PRIVATE END

Matrix *matrix_raw(int rows, int cols) {
  Matrix *matrix = malloc(sizeof(Matrix));
  matrix->rows = rows;
  matrix->cols = cols;
  matrix->data = malloc(rows * cols * sizeof(void *));
  return matrix;
}

Matrix *matrix_zeros(int rows, int cols) {
  Matrix *matrix = matrix_raw(rows, cols);
  void **data = matrix->data;
  for (int i = 0; i < rows * cols; ++i) {
    *((int *) (data[i] = malloc(sizeof(int)))) = 0;
  }
  return matrix;
}

Matrix *matrix_create(int rows, int cols, int value) {
  Matrix *matrix = matrix_raw(rows, cols);
  void **data = matrix->data;
  for (int i = 0; i < rows * cols; ++i) {
    *((int *) (data[i] = malloc(sizeof(int)))) = value;
  }
  return matrix;
}

Matrix *matrix_double_zeros(int rows, int cols) {
  Matrix *matrix = matrix_raw(rows, cols);
  void **data = matrix->data;
  for (int i = 0; i < rows * cols; ++i) {
    *((double *) (data[i] = malloc(sizeof(double)))) = 0;
  }
  return matrix;
}

Matrix *matrix_double_create(int rows, int cols, double value) {
  Matrix *matrix = matrix_raw(rows, cols);
  void **data = matrix->data;
  for (int i = 0; i < rows * cols; ++i) {
    *((double *) (data[i] = malloc(sizeof(double)))) = value;
  }
  return matrix;
}

Matrix *matrix_range(int from, int to) {
  Matrix *matrix = matrix_zeros((to - from) + 1, 1);
  for (int i = from; i <= to; ++i) {
    *((int *) matrix_element_by_index(matrix, i - from)) = i;
  }
  return matrix;
}

Matrix *matrix_from_file(char *file, bool transposed) {
  FILE *f = fopen(file, "r");
  assert(f != NULL);
  char num[11] = "";
  char c;
  int lc = 0, rows = 0, cols = 0;
  List *l_rows = list_empty();
  List *nums = list_empty();
  while ((c = fgetc(f)) != EOF) {
    assert(c == ',' || c == '\n' || (c >= '0' && c <= '9'));
    if (c == ',' || c == '\n') {
      assert(strlen(num) > 0);
      int val;
      sscanf(num, "%d", &val);
      list_push_int(nums, val);
      num[0] = '\0';
      ++cols;
      if (c == '\n') {
        ++rows;
        list_push(l_rows, nums);
        nums = list_empty();
        if (lc != 0) assert(cols == lc);
        lc = cols;
        cols = 0;
      }
    } else {
      char s[2] = {c, '\0'};
      strcat(num, s);
    }
  }
  list_delete(nums);
  fclose(f);
  Matrix *m = transposed ? matrix_zeros(l_rows->count, lc) : matrix_zeros(lc, l_rows->count);
  for (int r = 0; r < l_rows->count; ++r) {
    List *cv = list_get(l_rows, r);
    for (int c = 0; c < lc; ++c) {
      *(int *) matrix_element(m, transposed ? r : c, transposed ? c : r) = list_get_int(cv, c);
    }
    list_delete(cv);
  }
  list_scrap(l_rows);
  return m;
}

void matrix_to_file(Matrix *m, char *file) {
  FILE *f = fopen(file, "w");
  for (int r = 0; r < m->rows; ++r) {
    for (int c = 0; c < m->cols - 1; ++c) {
      fprintf(f, "%d,", *(int *) matrix_element(m, r, c));
    }
    fprintf(f, "%d\n", *(int *) matrix_element(m, r, m->cols - 1));
  }
  fclose(f);
}

void *matrix_element(Matrix *matrix, int row, int col) {
  if (row >= matrix->rows || col >= matrix->cols) {
    exit(EXIT_FAILURE);
  }
  return (matrix->data)[_matrix_index_for(matrix, row, col)];
}

void *matrix_element_by_index(Matrix *matrix, int index) {
  if (index >= matrix->rows * matrix->cols) {
    exit(EXIT_FAILURE);
  }
  return (matrix->data)[index];
}

void matrix_set(Matrix *matrix, int value) {
  int **data = (int **) matrix->data;
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    *data[i] = value;
  }
}

int matrix_prod(Matrix *matrix) {
  int p = 1;
  void **data = matrix->data;
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    p *= *(int *) data[i];
  }
  return p;
}

int matrix_sum(Matrix *matrix) {
  int p = 0;
  void **data = matrix->data;
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    p += *(int *) data[i];
  }
  return p;
}

double matrix_double_sum(Matrix *matrix) {
  double p = 0;
  void **data = matrix->data;
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    p += *(double *) data[i];
  }
  return p;
}

void matrix_double_set(Matrix *matrix, double value) {
  double **data = (double **) matrix->data;
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    *data[i] = value;
  }
}

Matrix *matrix_sub_indices(Matrix *o_matrix, int row_start, int row_end, int col_start, int col_end) {
  Matrix *n_matrix = matrix_raw(row_end - row_start, col_end - col_start);
  void **n_data = n_matrix->data, **o_data = o_matrix->data;
  for (int c = col_start, ci = 0; c < col_end; ++c, ++ci) {
    assert(c < o_matrix->cols);
    for (int r = row_start, ri = 0; r < row_end; ++r, ++ri) {
      assert(r < o_matrix->rows);
      n_data[_matrix_index_for(n_matrix, ri, ci)] = o_data[_matrix_index_for(o_matrix, r, c)];
    }
  }
  return n_matrix;
}

Matrix *matrix_sub_lists(Matrix *o_matrix, List *rows, List *cols) {
  Matrix *n_matrix = matrix_raw(rows->count, cols->count);
  void **n_data = n_matrix->data, **o_data = o_matrix->data;
  for (int ri = 0; ri < rows->count; ++ri) {
    for (int ci = 0, r = list_get_int(rows, ri); ci < cols->count; ++ci) {
      int c = list_get_int(cols, ci);
      assert(r < o_matrix->rows && c < o_matrix->cols);
      n_data[_matrix_index_for(n_matrix, ri, ci)] = o_data[_matrix_index_for(o_matrix, r, c)];
    }
  }
  return n_matrix;
}

Matrix *matrix_sub_list_index(Matrix *o_matrix, List *rows, int col_start, int col_end) {
  Matrix *n_matrix = matrix_raw(rows->count, col_end - col_start);
  void **n_data = n_matrix->data, **o_data = o_matrix->data;
  for (int ri = 0; ri < rows->count; ++ri) {
    for (int c = col_start, ci = 0, r = list_get_int(rows, ri); c < col_end; ++c, ++ci) {
      assert(r < o_matrix->rows && c < o_matrix->cols);
      n_data[_matrix_index_for(n_matrix, ri, ci)] = o_data[_matrix_index_for(o_matrix, r, c)];
    }
  }
  return n_matrix;
}

Matrix *matrix_sub_index_list(Matrix *o_matrix, int row_start, int row_end, List *cols) {
  Matrix *n_matrix = matrix_raw(row_end - row_start, cols->count);
  void **n_data = n_matrix->data, **o_data = o_matrix->data;
  for (int r = row_start, ri = 0; r < row_end; ++r, ++ri) {
    for (int ci = 0; ci < cols->count; ++ci) {
      int c = list_get_int(cols, ci);
      assert(r < o_matrix->rows && c < o_matrix->cols);
      n_data[_matrix_index_for(n_matrix, ri, ci)] = o_data[_matrix_index_for(o_matrix, r, c)];
    }
  }
  return n_matrix;
}

Matrix *matrix_sub_col(Matrix *matrix, int col) {
  return matrix_sub_indices(matrix, 0, matrix->rows, col, col + 1);
}

Matrix *matrix_sub_row(Matrix *matrix, int row) {
  return matrix_sub_indices(matrix, row, row + 1, 0, matrix->cols);
}

Matrix *matrix_sub_concat_rows(Matrix *matrix, Matrix *rows, bool parallel) {
  assert(matrix->cols == rows->cols);
  Matrix *m = matrix_raw(matrix->rows + rows->rows, matrix->cols);
  int **m1_data = (int **) matrix->data, **m2_data = (int **) rows->data,
      **m_data = (int **) m->data;
#pragma omp parallel for if (parallel) collapse(2)
  for (int r = 0; r < matrix->rows; ++r) {
    for (int c = 0; c < matrix->cols; ++c) {
      m_data[_matrix_index_for(m, r, c)] = m1_data[_matrix_index_for(matrix, r, c)];
    }
  }
#pragma omp parallel for if (parallel) collapse(2)
  for (int r = 0; r < rows->rows; ++r) {
    for (int c = 0; c < rows->cols; ++c) {
      m_data[_matrix_index_for(m, matrix->rows + r, c)] = m2_data[_matrix_index_for(rows, r, c)];
    }
  }
  return m;
}

List *matrix_find_by_value(Matrix *matrix, int value) {
  List *list = list_empty();
  for (int c = 0; c < matrix->cols; ++c) {
    for (int r = 0; r < matrix->rows; ++r) {
      if (*((int *) matrix_element(matrix, r, c)) == value) list_push_int(list, _matrix_index_for(matrix, r, c));
    }
  }
  return list;
}

List *matrix_double_find_by_value(Matrix *matrix, double value) {
  List *list = list_empty();
  for (int c = 0; c < matrix->cols; ++c) {
    for (int r = 0; r < matrix->rows; ++r) {
      if (*((double *) matrix_element(matrix, r, c)) == value) list_push_double(list, _matrix_index_for(matrix, r, c));
    }
  }
  return list;
}

List *matrix_to_list(Matrix *matrix) {
  List *l = list_empty();
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) list_push_int(l, *((int *) matrix_element_by_index(matrix, i)));
  return l;
}

List *matrix_double_to_list(Matrix *matrix) {
  List *l = list_empty();
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) list_push_double(l, *((double *) matrix_element_by_index(matrix, i)));
  return l;
}

Matrix *matrix_from_array(int rows, int cols, int arr[]) {
  Matrix *m = matrix_zeros(rows, cols);
  for (int i = 0; i < rows * cols; ++i) {
    *(int *) matrix_element_by_index(m, i) = arr[i];
  }
  return m;
}

Matrix *matrix_from_list(List *l) {
  Matrix *m = matrix_zeros(l->count, 1);
  for (int i = 0; i < l->count; ++i) {
    *(int *) matrix_element_by_index(m, i) = list_get_int(l, i);
  }
  return m;
}

void matrix_delete(Matrix *matrix) {
  for (int index = 0; index < matrix->rows * matrix->cols; ++index) {
    free((matrix->data)[index]);
  }
  matrix_scrap(matrix);
}

void matrix_scrap(Matrix *matrix) {
  free(matrix->data);
  matrix->data = NULL;
  free(matrix);
}

void *matrix_element_n_dim(Matrix *m, Matrix *ind, Matrix *dims) {
  assert(m->rows * m->cols == matrix_prod(dims));
  assert(ind->rows * ind->cols == dims->rows * dims->cols);
  int index = 0;
  int **id = (int **) ind->data, **dd = (int **) dims->data;
  for (int i = 0, m = 1; i < ind->rows * ind->cols; m *= *dd[i++]) {
    assert(*id[i] < *dd[i]);
    index += *id[i] * m;
  }
  return matrix_element_by_index(m, index);
}

void matrix_display(Matrix *matrix) {
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    printf("[%d]", *((int *) matrix_element_by_index(matrix, i)));
  }
  printf("\n");
}

void matrix_double_display(Matrix *matrix) {
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    printf("[%f]", *((double *) matrix_element_by_index(matrix, i)));
  }
  printf("\n");
}

Matrix *matrix_add(Matrix *a, Matrix *b) {
  assert(a->rows == b->rows);
  assert(a->cols == b->cols);
  Matrix *r = matrix_zeros(a->rows, a->cols);
  int **ad = (int **) a->data, **bd = (int **) b->data, **rd = (int **) r->data;
  for (int i = 0; i < a->rows * a->cols; ++i) {
    *rd[i] = *ad[i] + *bd[i];
  }
  return r;
}

void matrix_add_in(Matrix *a, Matrix *b) {
  assert(a->rows == b->rows);
  assert(a->cols == b->cols);
  int **ad = (int **) a->data, **bd = (int **) b->data;
  for (int i = 0; i < a->rows * a->cols; ++i) {
    *ad[i] = *ad[i] + *bd[i];
  }
}

Matrix *matrix_double_add(Matrix *a, Matrix *b) {
  assert(a->rows == b->rows);
  assert(a->cols == b->cols);
  Matrix *r = matrix_double_zeros(a->rows, a->cols);
  double **ad = (double **) a->data, **bd = (double **) b->data, **rd = (double **) r->data;
  for (int i = 0; i < a->rows * a->cols; ++i) {
    *rd[i] = *ad[i] + *bd[i];
  }
  return r;
}

Matrix *matrix_add_int_double(Matrix *a, Matrix *b) {
  assert(a->rows == b->rows);
  assert(a->cols == b->cols);
  Matrix *r = matrix_double_zeros(a->rows, a->cols);
  int **ad = (int **) a->data;
  double **bd = (double **) b->data, **rd = (double **) r->data;
  for (int i = 0; i < a->rows * a->cols; ++i) {
    *rd[i] = *ad[i] + *bd[i];
  }
  return r;
}

Matrix *matrix_sum_n_cols(Matrix *matrix, int cols) {
  assert((matrix->rows * matrix->cols) % cols == 0);
  Matrix *m = matrix_zeros((matrix->rows * matrix->cols) / cols, 1);
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    *(int *) matrix_element_by_index(m, i % m->rows) += *(int *) matrix_element_by_index(matrix, i);
  }
  return m;
}

Matrix *matrix_double_sum_n_cols(Matrix *matrix, int cols) {
  assert((matrix->rows * matrix->cols) % cols == 0);
  Matrix *m = matrix_double_zeros((matrix->rows * matrix->cols) / cols, 1);
  for (int i = 0; i < matrix->rows * matrix->cols; ++i) {
    *(double *) matrix_element_by_index(m, i % m->rows) += *(double *) matrix_element_by_index(matrix, i);
  }
  return m;
}

Matrix *matrix_double_subtract(Matrix *m1, Matrix *m2) {
  Matrix *m = matrix_double_zeros(m1->rows, m1->cols);
  for (int i = 0; i < m1->rows * m1->cols; ++i) {
    *(double *) matrix_element_by_index(m, i) = *(double *) matrix_element_by_index(m1, i) - *(double *) matrix_element_by_index(m2, i);
  }
  return m;
}

Matrix *matrix_lgamma(Matrix *m) {
  Matrix *nm = matrix_double_zeros(m->rows, m->cols);
  for (int i = 0; i < m->rows * m->cols; ++i) {
    *(double *) matrix_element_by_index(nm, i) = lgamma(*(double *) matrix_element_by_index(m, i));
  }
  return nm;
}

Matrix *matrix_create_sz(Matrix *matrix) {
  Matrix *r = matrix_zeros(matrix->rows, 1);
  for (int i = 0; i < matrix->rows; ++i) {
    int max = 2;
    for (int temp = 1, j = 0; j < matrix->cols; ++j) {
      temp = *(int *) matrix_element(matrix, i, j);
      max = max < temp ? temp : max;
    }
    *(int *) matrix_element_by_index(r, i) = max;
  }
  return r;
}

void matrix_double_mk_stochastic(Matrix *m, Matrix *ns) {
  int dim = *(int *) matrix_element_by_index(ns, (ns->rows * ns->cols) - 1);
  int index_count = m->rows * m->cols;
  assert(index_count % dim == 0);
  Matrix *div = matrix_double_sum_n_cols(m, dim);
  for (int i = 0; i < index_count; ++i) {
    *(double *) matrix_element_by_index(m, i) = (*(double *) matrix_element_by_index(m, i)) / (*(double *) matrix_element_by_index(div, i % (index_count / dim)));
  }
  matrix_delete(div);
}
