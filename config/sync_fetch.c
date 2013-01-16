#include <stdlib.h>
int main () {
  int o;
  __sync_fetch_and_add(&o, 1);
  __sync_sub_and_fetch(&o, 1);
  exit (0);
}
