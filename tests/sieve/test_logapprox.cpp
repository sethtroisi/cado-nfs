#include "cado.h"
#include <string.h>
#include <iostream>
#include "logapprox.hpp"

void display_logapprox(piecewise_linear_function & F)
{
    std::cout << "Total " << F.equations.size() << " pieces\n";
    for( ; !F.equations.empty() ; ) {
        double r0 = F.endpoints.front();
        F.endpoints.pop_front();
        double r1 = F.endpoints.front();
        std::pair<double,double> uv = F.equations.front();
        F.equations.pop_front();
        std::cout << "[" << r0 << ", " << r1 << "]: " << uv.first << " + x * " << uv.second << "\n";
    }
}

int test_from_bug21684(bool display)
{
    cxx_double_poly f(3);
    for(double x : {-6.3406659802246472e+28, 6.4695148868405632e+28, 9.5457310557271272e+27}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.3300700859809263);
    piecewise_linear_function F = A.logapprox(-2048,2048);
    ASSERT_ALWAYS(F.equations.size() == 21);
    if (display) display_logapprox(F);
    return 0;
}

int original_test(bool display)
{
    cxx_double_poly f(6);
    for(double x : {-39769265437440000., -302859053854976., 5377439145928., -1684314626., -5377481., 970., 1.}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.1);
    piecewise_linear_function F = A.logapprox(-2048,2048);
    ASSERT_ALWAYS(F.equations.size() == 34);
    if (display) display_logapprox(F);
    return 0;
}

int main(int argc, char * argv[])
{
    bool display = false;
    if (argc >= 2 && strcmp(argv[1], "--display") == 0)
        display=true;
    original_test(display);
    test_from_bug21684(display);
}

