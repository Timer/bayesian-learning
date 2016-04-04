#ifndef LIST_H
#define LIST_H

typedef struct {
  int size, count;
  void **arr;
} List;

List *list_empty();
void list_enlarge(List *list);
void list_grow(List *list, int size);
int list_push(List *list, void *data);
int list_push_int(List *list, int value);
int list_push_double(List *list, double value);
void *list_set(List *list, int index, void *data);
void *list_get(List *list, int index);
int list_get_int(List *list, int index);
double list_get_double(List *list, int index);
void *list_remove(List *list, int index);
List *list_slice(List *list, int start, int end);
List *list_difference_type_int(List *primary, List *secondary);
void list_delete(List *list);
void list_scrap(List *list);

#endif
