#include "bnet.h"
#include <assert.h>
#include <stdlib.h>
#include "matrix.h"

void cpd_delete(CPD *cpd) {
  matrix_delete(cpd->sizes);
  matrix_delete(cpd->dirichlet);
  if (cpd->cpt != NULL) matrix_delete(cpd->cpt);
  free(cpd);
}

List *adjacency_matrix_parents(Matrix *adj_mat, int col) {
  List *l = list_empty();
  Matrix *sub = matrix_sub_indices(adj_mat, 0, adj_mat->rows, col, col + 1);
  for (int i = 0; i < adj_mat->rows; ++i) {
    int val = *(int *) matrix_element_by_index(sub, i);
    assert(val == 0 || val == 1);
    if (val) list_push_int(l, i);
  }
  matrix_scrap(sub);
  return l;
}
