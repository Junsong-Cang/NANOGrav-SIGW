#include <stdio.h>
#include <gsl/gsl_sf_hyperg.h>

int main() {
    double a = 1.5;  // Parameter a
    double b = 2.5;  // Parameter b
    double x = 0.8;  // Argument x

    double result = gsl_sf_hyperg_1F1(a, b, x);
    printf("Result: %g\n", result);

    return 0;
}
