#include "./misc/print_generic_sgf.hpp"

#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    // Seawater example from p. 23 in Michalski2016b.
    FrequencyDomain fd(10e6);
    cmplx eps_r = cmplx_permittivity(fd, 81.0, 4.0);

    arma::vec x_vals = arma::linspace(0, 10, 70) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-1, 7, 70) * fd.lambda_0;

    print_generic_sgf(fd, eps_r, x_vals, z_vals);

    return 0;
}
