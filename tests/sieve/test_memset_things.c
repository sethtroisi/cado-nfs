#include "cado.h"
#include "las-norms.h"
#include "utils.h"

int main()
{
    verbose_output_init(3);
    verbose_output_add(0, stdout, 2);
    tune_las_memset();
    return 0;
}

