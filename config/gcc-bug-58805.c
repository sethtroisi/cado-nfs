/* 
  This program is a test case for GCC bug 58805, see
  http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58805
  This bug affects the cofactorization routines.
  If the gcc version is unaffected, the compiled executable returns with
  exit status 0, otherwise it aborts.
*/

static inline void bar(unsigned long *r)
{
  unsigned long t;
  __asm__ (
    "movq $42, %[t]\n\t"
    "movq %[t], %[r]\n\t"
    : [t] "=&r" (t), [r] "=r" (*r)
  );
}

void foo(int n, unsigned long *x, unsigned long *y)
{
  if (n == 0)
    bar(x);
  else
    bar(y);
}

int n = 1;
unsigned long x, y;

int main (void)
{
  foo(n, &x, &y);
  if (y != 42)
    __builtin_abort ();
  return 0;
}
