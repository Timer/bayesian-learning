#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "bnet.h"
#include "list.h"
#include "matrix.h"
#include "rand.h"

#define MPI 0
#define VERBOSE 0
#define SAVE_NETWORKS 0

#if MPI
#include <mpi.h>
#endif

double dirichlet_score_family(Matrix *counts, CPD *cpd) {
  Matrix *ns = cpd->sizes, *prior = cpd->dirichlet;
  Matrix *ns_self = matrix_sub_indices(ns, ns->rows - 1, ns->rows, 0, ns->cols);
  Matrix *pnc = matrix_add_int_double(counts, prior);
  Matrix *gamma_pnc = matrix_lgamma(pnc), *gamma_prior = matrix_lgamma(prior);
  matrix_delete(pnc);
  Matrix *lu_mat = matrix_double_subtract(gamma_pnc, gamma_prior);
  matrix_delete(gamma_pnc);
  matrix_delete(gamma_prior);
  Matrix *LU = matrix_double_sum_n_cols(lu_mat, *(int *) matrix_element_by_index(ns_self, 0));
  matrix_delete(lu_mat);
  Matrix *alpha_ij = matrix_double_sum_n_cols(prior, *(int *) matrix_element_by_index(ns_self, 0));
  Matrix *N_ij = matrix_sum_n_cols(counts, *(int *) matrix_element_by_index(ns_self, 0));
  matrix_scrap(ns_self);
  Matrix *gamma_alpha = matrix_lgamma(alpha_ij);
  Matrix *alpha_N = matrix_add_int_double(N_ij, alpha_ij);
  matrix_delete(N_ij);
  matrix_delete(alpha_ij);
  Matrix *gamma_alpha_N = matrix_lgamma(alpha_N);
  matrix_delete(alpha_N);
  Matrix *LV = matrix_double_subtract(gamma_alpha, gamma_alpha_N);
  matrix_delete(gamma_alpha);
  matrix_delete(gamma_alpha_N);
  Matrix *LU_LV = matrix_double_add(LU, LV);
  matrix_delete(LU);
  matrix_delete(LV);
  double score = matrix_double_sum(LU_LV);
  matrix_delete(LU_LV);
  return score;
}

int count_index(Matrix *sz, Matrix *sample_data, int col) {
  Matrix *mat_col = matrix_sub_col(sample_data, col);
  int index = 0;
  int **id = (int **) mat_col->data, **dd = (int **) sz->data;
  for (int i = 0, m = 1; i < mat_col->rows * mat_col->cols; m *= *dd[i++]) {
    assert((*id[i]) - 1 < *dd[i]);
    index += ((*id[i]) - 1) * m;
  }
  matrix_scrap(mat_col);
  return index;
}

Matrix *compute_counts(Matrix *data, Matrix *sz) {
  assert(sz->rows == data->rows);
  Matrix *count = matrix_zeros(matrix_prod(sz), 1);
  for (int i = 0; i < data->cols; ++i) {
    *((int *) matrix_element_by_index(count, count_index(sz, data, i))) += 1;
  }
  return count;
}

double log_marg_prob_node(CPD *cpd, Matrix *self_ev, Matrix *pev) {
  assert(self_ev->rows == 1);
  assert(cpd->sizes->cols == 1);
  assert(pev->cols == self_ev->cols);
  Matrix *data = matrix_sub_concat_rows(pev, self_ev, false);
  Matrix *counts = compute_counts(data, cpd->sizes);
  matrix_scrap(data);
  double score = dirichlet_score_family(counts, cpd);
  matrix_delete(counts);
  return score;
}

Matrix *prob_node(CPD *cpd, Matrix *self_ev, Matrix *pev) {
  Matrix *sample_data = matrix_sub_concat_rows(pev, self_ev, false);
  Matrix *prob = matrix_double_zeros(sample_data->rows, sample_data->cols);
  for (int i = 0; i < sample_data->cols; ++i) {
    Matrix *mat_col = matrix_sub_col(sample_data, i);
    int index = 0;
    int **id = (int **) mat_col->data, **dd = (int **) cpd->sizes->data;
    for (int j = 0, m = 1; j < mat_col->rows * mat_col->cols; m *= *dd[j++]) {
      assert((*id[j]) - 1 < *dd[j]);
      index += ((*id[j]) - 1) * m;
    }
    matrix_scrap(mat_col);
    *(double *) matrix_element_by_index(prob, i) = *(double *) matrix_element_by_index(cpd->cpt, index);
  }
  matrix_scrap(sample_data);
  return prob;
}

double log_prob_node(CPD *cpd, Matrix *self_ev, Matrix *pev) {
  double score = 0;
  Matrix *p = prob_node(cpd, self_ev, pev);
  double **data = (double **) p->data;
  for (int i = 0; i < p->rows * p->cols; ++i, ++data) {
    double d = **data;
    score += d <= 0 ? DBL_MIN : log(d);
  }
  matrix_delete(p);
  return score;
}

CPD *tabular_CPD(Matrix *dag, Matrix *ns, int self) {
  CPD *cpd = malloc(sizeof(CPD));
  List *ps = adjacency_matrix_parents(dag, self);
  list_push_int(ps, self);
  Matrix *fam_sz = matrix_zeros(ps->count, 1);
  for (int i = 0; i < ps->count; ++i) {
    *(int *) matrix_element_by_index(fam_sz, i) = *(int *) matrix_element_by_index(ns, list_get_int(ps, i));
  }
  cpd->sizes = fam_sz;
  Matrix *calc = matrix_sub_indices(fam_sz, 0, ps->count - 1, 0, 1);
  int psz = matrix_prod(calc);
  list_delete(ps);
  matrix_scrap(calc);
  cpd->dirichlet = matrix_double_create(matrix_prod(fam_sz), 1, (1.0 / psz) * (1.0 / *(int *) matrix_element_by_index(ns, self)));
  cpd->cpt = NULL;
  return cpd;
}

double score_family(int j, List *ps, Matrix *ns, List *discrete, Matrix *data, char *scoring_fn) {
  Matrix *dag = matrix_zeros(data->rows, data->rows);
  if (ps->count > 0) {
    Matrix *dag_sub = matrix_sub_list_index(dag, ps, j, j + 1);
    matrix_set(dag_sub, 1);
    matrix_scrap(dag_sub);
    //TODO: sort `ps` here.
  }
  Matrix *data_sub_1 = matrix_sub_indices(data, j, j + 1, 0, data->cols),
         *data_sub_2 = matrix_sub_list_index(data, ps, 0, data->cols);
  CPD *cpd = tabular_CPD(dag, ns, j);
  double score;
  if (!strcmp(scoring_fn, "bayesian")) {
    score = log_marg_prob_node(cpd, data_sub_1, data_sub_2);
  } else if (!strcmp(scoring_fn, "bic")) {
    List *fam = list_slice(ps, 0, ps->count);
    int a_index = list_push_int(fam, j);
    Matrix *data_sub_3 = matrix_sub_list_index(data, fam, 0, data->cols);
    Matrix *counts = compute_counts(data_sub_3, cpd->sizes);
    matrix_scrap(data_sub_3);
    cpd->cpt = matrix_add_int_double(counts, cpd->dirichlet);
    matrix_delete(counts);
    matrix_double_mk_stochastic(cpd->cpt, ns);
    double L = log_prob_node(cpd, data_sub_1, data_sub_2);
    Matrix *sz = cpd->sizes;
    int *last = (int *) sz->data[sz->rows * sz->cols - 1];
    --*last;
    score = L - 0.5 * matrix_prod(sz) * log(data->cols);
    ++*last;
    free(list_remove(fam, a_index));
    list_scrap(fam);
  } else {
    assert(1 == 2);
  }

  cpd_delete(cpd);
  matrix_scrap(data_sub_1);
  matrix_scrap(data_sub_2);
  matrix_delete(dag);
  return score;
}

Matrix *learn_struct_K2(Matrix *data, Matrix *ns, List *order, char *scoring_fn, int max_parents) {
  assert(order->count == data->rows);
  int n = data->rows;
  int max_fan_in = max_parents == 0 ? n : max_parents;
  List *discrete = list_empty();
  for (int i = 0; i < n; ++i) list_push_int(discrete, i);

  Matrix *dag = matrix_zeros(n, n);
  int parent_order = 0;
  for (int i = 0; i < n; ++i) {
    List *ps = list_empty();
    int j = list_get_int(order, i);
    double score = score_family(j, ps, ns, discrete, data, scoring_fn);
#if VERBOSE
    printf("\nnode %d, empty score %6.4f\n", j, score);
#endif
    for (; ps->count <= max_fan_in;) {
      List *order_sub = list_slice(order, 0, i);
      List *pps = list_difference_type_int(order_sub, ps);
      list_scrap(order_sub);
      int nps = pps->count;
      Matrix *pscore = matrix_double_zeros(1, nps);
      for (int pi = 0; pi < nps; ++pi) {
        int p = list_get_int(pps, pi);
        int n_index = list_push_int(ps, p);
        *((double *) matrix_element_by_index(pscore, pi)) = score_family(j, ps, ns, discrete, data, scoring_fn);
#if VERBOSE
        printf("considering adding %d to %d, score %6.4f\n", p, j, *((double *) matrix_element_by_index(pscore, pi)));
#endif
        free(list_remove(ps, n_index));
      }
      double best_pscore = -DBL_MAX;
      int best_p = -1;
      for (int i = 0; i < nps; ++i) {
        double d = *(double *) matrix_element_by_index(pscore, i);
        if (d > best_pscore) {
          best_pscore = d;
          best_p = i;
        }
      }
      matrix_delete(pscore);
      if (best_p == -1) {
        list_scrap(pps);
        break;
      }
      best_p = list_get_int(pps, best_p);
      list_scrap(pps);
      if (best_pscore > score) {
        score = best_pscore;
        list_push_int(ps, best_p);
#if VERBOSE
        printf("* adding %d to %d, score %6.4f\n", best_p, j, best_pscore);
#endif
      } else {
        break;
      }
    }
    if (ps->count > 0) {
      Matrix *dag_sub = matrix_sub_list_index(dag, ps, j, j + 1);
      matrix_set(dag_sub, ++parent_order);
      matrix_scrap(dag_sub);
    }
    list_delete(ps);
  }
  list_delete(discrete);
  return dag;
}

#if MPI
#define MPI_TAG_MATRIX_R 1
#define MPI_TAG_MATRIX_C 2
#define MPI_TAG_MATRIX_D 3

void MPI_Matrix_Send(int to_index, Matrix *matrix) {
  int rows = matrix->rows, cols = matrix->cols;
  int **data = (int **) matrix->data;
  int *extracted = malloc(rows * cols * sizeof(int));
  int *extracted_start = extracted;
  for (int i = 0; i < rows * cols; ++i, ++extracted, ++data) {
    *extracted = **data;
  }
  MPI_Send(&rows, 1, MPI_INT, to_index, MPI_TAG_MATRIX_R, MPI_COMM_WORLD);
  MPI_Send(&cols, 1, MPI_INT, to_index, MPI_TAG_MATRIX_C, MPI_COMM_WORLD);
  MPI_Send(extracted_start, rows * cols, MPI_INT, to_index, MPI_TAG_MATRIX_D, MPI_COMM_WORLD);
  free(extracted_start);
}

Matrix *MPI_Matrix_Recv(int from_index) {
  int rows, cols;
  MPI_Recv(&rows, 1, MPI_INT, from_index, MPI_TAG_MATRIX_R, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&cols, 1, MPI_INT, from_index, MPI_TAG_MATRIX_C, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  Matrix *matrix = matrix_zeros(rows, cols);
  int **data = (int **) matrix->data;
  int *values = malloc(rows * cols * sizeof(int));
  MPI_Recv(values, rows * cols, MPI_INT, from_index, MPI_TAG_MATRIX_D, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  int *values_start = values;
  for (int i = 0; i < rows * cols; ++i, ++values, ++data) {
    **data = *values;
  }
  free(values_start);
  return matrix;
}
#endif

int exec(int forkIndex, int forkSize, bool data_transposed, char *f_data, int topologies, char *f_output, char *scoring_fn, int max_parents) {
  Matrix *data = matrix_from_file(f_data, data_transposed), *sz = matrix_create_sz(data);
#if MPI
  assert(forkIndex > -1);
  assert(forkSize > 0);
  int top_d = topologies / forkSize, top_r = topologies % forkSize;
  if (forkIndex < top_r) ++top_d;
  topologies = top_d;
#endif
  Matrix *orders = matrix_zeros(data->rows * topologies, data->rows);
#if SAVE_NETWORKS
  Matrix *networks = matrix_zeros(data->rows * topologies, data->rows * data->rows);
#endif

#pragma omp parallel for
  for (int r = 0; r < orders->rows; ++r) {
    int start = r / topologies;
    int *arr = malloc(orders->cols * sizeof(int));
    arr[0] = start;
    for (int i = 1; i < orders->cols; ++i) {
      arr[i] = i == start ? 0 : i;
    }
    shuffle_int(orders->cols - 1, arr + 1);
    for (int c = 0; c < orders->cols; ++c) {
      *(int *) matrix_element(orders, r, c) = arr[c];
    }
    free(arr);
  }

  Matrix *consensus_network = matrix_zeros(data->rows, data->rows);
  int cn_n_elements = consensus_network->rows * consensus_network->cols;

#pragma omp parallel for
  for (int o = 0; o < orders->rows; ++o) {
    Matrix *m_order = matrix_sub_indices(orders, o, o + 1, 0, orders->cols);
    List *order = matrix_to_list(m_order);
    Matrix *bnet = learn_struct_K2(data, sz, order, scoring_fn, max_parents);
    assert(consensus_network->rows == bnet->rows);
    assert(consensus_network->cols == bnet->cols);

#pragma omp critical
    for (int i = 0; i < cn_n_elements; ++i) {
      *(int *) matrix_element_by_index(consensus_network, i) += *(int *) matrix_element_by_index(bnet, i) ? 1 : 0;
    }

#if SAVE_NETWORKS
    for (int i = 0; i < cn_n_elements; ++i) {
      *(int *) matrix_element_by_index(networks, i + cn_n_elements * o) = *(int *) matrix_element_by_index(bnet, i);
    }
#endif

    matrix_delete(bnet);
    list_delete(order);
    matrix_scrap(m_order);
  }

  matrix_delete(sz);
  matrix_delete(data);

#if MPI
  //TODO: merge and write topologies
  if (forkIndex == 0) {
    for (int i = 1; i < forkSize; ++i) {
      Matrix *merge = MPI_Matrix_Recv(i);
      matrix_add_in(consensus_network, merge);
      matrix_delete(merge);
    }
#if SAVE_NETWORKS
    for (int i = 1; i < forkSize; ++i) {
      Matrix *merge = MPI_Matrix_Recv(i), *old = orders;
      orders = matrix_sub_concat_rows(orders, merge, true);
      matrix_scrap(old);
      matrix_scrap(merge);
    }
    for (int i = 1; i < forkSize; ++i) {
      Matrix *merge = MPI_Matrix_Recv(i), *old = networks;
      networks = matrix_sub_concat_rows(networks, merge, true);
      matrix_scrap(old);
      matrix_scrap(merge);
    }
#endif
  } else {
    MPI_Matrix_Send(0, consensus_network);
#if SAVE_NETWORKS
    MPI_Matrix_Send(0, orders);
    MPI_Matrix_Send(0, networks);
#endif
  }
#endif
  if (forkIndex == 0) {
    matrix_to_file(consensus_network, f_output);
#if SAVE_NETWORKS
    matrix_to_file(networks, "networks.csv");
    matrix_to_file(orders, "topologies.csv");
#endif
  }
#if SAVE_NETWORKS
  matrix_delete(networks);
#endif
  matrix_delete(orders);
  matrix_delete(consensus_network);
  return 0;
}

int main(int argc, char **argv) {
  int forkIndex = 0, forkSize = 1;
#if MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &forkIndex);
  MPI_Comm_size(MPI_COMM_WORLD, &forkSize);
#endif

  srand(time(NULL) ^ forkIndex);
  int threads = 1, topologies = 1, max_parents = 0;
  bool data_transposed = false;
  char *data = NULL, *output = "consensus.csv";
  char *scoring_fn = "bayesian";
  int c;
  while ((c = getopt(argc, argv, "Thp:d:t:o:s:")) != -1) {
    switch (c) {
    case 'T': {
      data_transposed = true;
      break;
    }
    case 'p': {
      threads = atoi(optarg);
      assert(threads > 0);
      assert(threads <= omp_get_num_procs());
      break;
    }
    case 'm': {
      max_parents = atoi(optarg);
      assert(max_parents >= 0);
      break;
    }
    case 'd': {
      data = optarg;
      break;
    }
    case 't': {
      topologies = atoi(optarg);
      break;
    }
    case 'o': {
      output = optarg;
      break;
    }
    case 's': {
      scoring_fn = optarg;
      break;
    }
    case 'h':
    default: {
      puts(": -p <num_threads> -d <data file> -t <topologies per gene> -o <output file> -m <max parents>");
      puts("~ -T (reads matrix transposed)");
      return 1;
    }
    }
  }
  if (data == NULL) {
    puts("You must send a data file using -d <file name>.");
    return 1;
  }
  omp_set_num_threads(threads);
  int status = exec(forkIndex, forkSize, data_transposed, data, topologies, output, scoring_fn, max_parents);
#if MPI
  MPI_Finalize();
#endif
  return status;
}
