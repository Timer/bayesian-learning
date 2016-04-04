#include "list.h"
#include <assert.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

List *list_empty() {
  List *list = malloc(sizeof(List));
  list->size = list->count = 0;
  list->arr = NULL;
  return list;
}

void list_enlarge(List *list) {
  if (list->arr) {
    if (list->count == list->size) {
      int prev = list->size;
      list->arr = realloc(list->arr, (list->size *= 2) * sizeof(void *));
      for (int i = prev; i < list->size; ++i) list_set(list, i, NULL);
    }
  } else {
    list->size = 2;
    list->arr = malloc(sizeof(void *) * list->size);
    for (int i = 0; i < list->size; ++i) list_set(list, i, NULL);
  }
}

void list_grow(List *list, int size) {
  while (list->count < size) list_push(list, NULL);
}

int list_push(List *list, void *data) {
  list_enlarge(list);
  (list->arr)[list->count] = data;
  return list->count++;
}

int list_push_int(List *list, int value) {
  list_enlarge(list);
  int *v = malloc(sizeof(int));
  *v = value;
  (list->arr)[list->count] = v;
  return list->count++;
}

int list_push_double(List *list, double value) {
  list_enlarge(list);
  int *v = malloc(sizeof(double));
  *v = value;
  (list->arr)[list->count] = v;
  return list->count++;
}

void *list_set(List *list, int index, void *data) {
  assert(index < list->size);
  void *old = (list->arr)[index];
  (list->arr)[index] = data;
  return old;
}

void *list_get(List *list, int index) {
  assert(index < list->count);
  return (list->arr)[index];
}

int list_get_int(List *list, int index) {
  assert(index < list->count);
  return *((int *) ((list->arr)[index]));
}

double list_get_double(List *list, int index) {
  assert(index < list->count);
  return *((double *) ((list->arr)[index]));
}

void *list_remove(List *list, int index) {
  assert(index < list->count--);
  void *ptr = (list->arr)[index];
  if (list->count - index > 0) {
    memmove(list->arr + index, list->arr + index + 1, (list->count - index) * sizeof(void *));
  }
  list->arr[list->count] = NULL;
  return ptr;
}

List *list_slice(List *ol, int start, int end) {
  List *sl = list_empty();
  for (int index = start; index < end; ++index) {
    list_push(sl, (ol->arr)[index]);
  }
  return sl;
}

List *list_difference_type_int(List *primary, List *secondary) {
  List *list = list_slice(primary, 0, primary->count);
  if (list->count == 0 || secondary->count == 0) return list;
  for (int i = 0; i < secondary->count; ++i) {
    int val = *((int *) (secondary->arr)[i]);
    for (int i2 = 0; i2 < list->count; ++i2) {
      if (*((int *) (list->arr)[i2]) == val) {
        list_remove(list, i2);
        break;
      }
    }
  }
  return list;
}

void list_delete(List *list) {
  for (int i = 0; i < list->size; ++i) {
    if (list->arr[i] == NULL) continue;
    free(list->arr[i]);
  }
  list_scrap(list);
}

void list_scrap(List *list) {
  free(list->arr);
  list->arr = NULL;
  free(list);
}
