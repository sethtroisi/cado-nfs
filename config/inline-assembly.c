
#include <stdlib.h>

static void 
addmul1_nc_1(unsigned long *z, const unsigned long *x, const unsigned long c)
{
  __asm__ volatile(
   "    movq    %2, %%rax\n"
   "    mulq    %[mult]\n"
   "    addq    %%rax, %0\n"
   "    adcq    $0, %%rdx\n"
   "    movq    %%rdx, %%rcx\n"
   "    addq    %%rcx, %1\n"
  : "+rm" (z[0]), "+rm" (z[1])
  : "rm" (x[0]), [mult] "r" (c)
  : "%rax", "%rcx", "%rdx");
}

int main()
{
    unsigned long z[2]={0,0};
    unsigned long x[1]={0};
    unsigned long c = 0;

    x[0] = rand(); c=rand(); addmul1_nc_1(z, x, c);
    x[0] = rand(); c=rand(); addmul1_nc_1(z, x, c);
    x[0] = rand(); c=rand(); addmul1_nc_1(z, x, c);
    x[0] = rand(); c=rand(); addmul1_nc_1(z, x, c);

    return z[0];
}

