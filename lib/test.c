#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "matrix.h"

int main(int argc, char **argv) {
  // --- MATRIX START
  puts("Creating test matrix ...");
  Matrix *m = matrix_zeros(3, 4);

  puts("Validating matrix set ...");
  matrix_set(m, 10);
  for (int i = 0; i < m->rows * m->cols; ++i) {
    assert(*(int *) matrix_element_by_index(m, i) == 10);
  }

  puts("Validating element and index correspondence ...");
  int test_row = 1, test_col = 2, test_index = 7;
  *(int *) matrix_element(m, test_row, test_col) *= 2;
  assert(*(int *) matrix_element_by_index(m, test_index) == 20);

  puts("Testing matrix find ...");
  List *l = matrix_find_by_value(m, 20);
  assert(l->count == 1);
  assert(list_get_int(l, 0) == test_index);
  list_delete(l);

  puts("Validating proper generation of sub matrices ...");
  Matrix *sm = matrix_sub_indices(m, 1, 2, 1, 3);
  assert(sm->rows == 1 && sm->cols == 2);
  matrix_set(sm, 5);
  matrix_scrap(sm);

  l = matrix_find_by_value(m, 5);
  assert(l->count == 2);
  assert(list_get_int(l, 0) == 4);
  assert(list_get_int(l, 1) == 7);
  list_delete(l);

  List *rows = list_empty(), *cols = list_empty();
  list_push_int(rows, 1);
  list_push_int(rows, 2);
  list_push_int(cols, 3);
  sm = matrix_sub_lists(m, rows, cols);
  assert(sm->rows == 2 && sm->cols == 1);
  matrix_set(sm, 6);
  matrix_scrap(sm);

  l = matrix_find_by_value(m, 6);
  assert(l->count == 2);
  assert(list_get_int(l, 0) == 10);
  assert(list_get_int(l, 1) == 11);
  list_delete(l);

  sm = matrix_sub_index_list(m, 1, 3, cols);
  assert(sm->rows == 2 && sm->cols == 1);
  matrix_set(sm, 7);
  matrix_scrap(sm);

  l = matrix_find_by_value(m, 7);
  assert(l->count == 2);
  assert(list_get_int(l, 0) == 10);
  assert(list_get_int(l, 1) == 11);
  list_delete(l);

  sm = matrix_sub_list_index(m, rows, 3, 4);
  assert(sm->rows == 2 && sm->cols == 1);
  matrix_set(sm, 18);
  matrix_scrap(sm);

  l = matrix_find_by_value(m, 18);
  assert(l->count == 2);
  assert(list_get_int(l, 0) == 10);
  assert(list_get_int(l, 1) == 11);
  list_delete(l);

  puts("Testing matrix_prod ...");
  matrix_set(m, 2);
  assert(matrix_prod(m) == 4096);

  matrix_delete(m);

  puts("Testing n-dimensional functionality ...");
  m = matrix_range(1, 8);
  int dims[] = {2, 2, 2};
  Matrix *m_ind = matrix_zeros(3, 1), *m_dims = matrix_from_array(3, 1, dims);
  //0, 0, 0
  assert(*(int *) matrix_element_n_dim(m, m_ind, m_dims) == 1);
  *(int *) matrix_element_by_index(m_ind, 0) = 1;
  //1, 0, 0
  assert(*(int *) matrix_element_n_dim(m, m_ind, m_dims) == 2);
  *(int *) matrix_element_by_index(m_ind, 1) = 1;
  //1, 1, 0
  assert(*(int *) matrix_element_n_dim(m, m_ind, m_dims) == 4);
  *(int *) matrix_element_by_index(m_ind, 1) = 0;
  *(int *) matrix_element_by_index(m_ind, 2) = 1;
  //1, 0, 1
  assert(*(int *) matrix_element_n_dim(m, m_ind, m_dims) == 6);
  *(int *) matrix_element_by_index(m_ind, 1) = 1;
  //1, 1, 1
  assert(*(int *) matrix_element_n_dim(m, m_ind, m_dims) == 8);

  matrix_delete(m_ind);
  matrix_delete(m_dims);
  matrix_delete(m);

  //TODO: test remaining matrix funcs
  //TODO: test count index / compute count functs
  // --- MATRIX END

  puts("All tests passed!");
  return 0;
}
