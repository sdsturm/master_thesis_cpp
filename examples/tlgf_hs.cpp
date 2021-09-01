#include <mthesis/mthesis.hpp>

#include <gsl/gsl_const_mksa.h>
#include <armadillo>

#include <iostream>

using namespace mthesis;

int main()
{
    auto fd = FrequencyDomain(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    auto ground = Medium(fd, 3, 1);
    auto hs = HalfSpace(fd, ground);

    real k_rho = 0.0;
    EmMode type = EmMode::TM;

    real z_ = fd.lambda_0;

    arma::vec z = arma::linspace(-1, 1, 1e3) * fd.lambda_0;
    arma::cx_vec V(z.n_elem), I(z.n_elem), Z(z.n_elem);

    for (arma::uword n = 0; n < z.n_elem; n++)
    {
        V(n) = V_i(hs, k_rho, z(n), z_, type);
        I(n) = I_i(hs, k_rho, z(n), z_, type);
        Z(n) = V(n) / I(n);
    }

    return 0;
}
