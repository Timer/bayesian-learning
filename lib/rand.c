#include <assert.h>
#include <stdlib.h>

unsigned int rand_inclusive(unsigned int min, unsigned int max) {
  assert(max >= min);
  unsigned int
      delta = 1 + max - min,
      b = RAND_MAX / delta, c = b * delta;
  int r;
  do {
    r = rand();
  } while (r >= c);
  return min + r / b;
}

void shuffle_int(int c, int *a) {
  int b, d;
  while (c) b = rand_inclusive(0, --c), d = a[c], a[c] = a[b], a[b] = d;
}
