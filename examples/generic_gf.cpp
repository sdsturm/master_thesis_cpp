#include <mthesis.hpp>

#include <gsl/gsl_const_mksa.h>
#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    // Problem specification.
    auto fd = FrequencyDomain(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    auto eps_r = cmplx_permittivity(fd, 3, 10e-3);
    // real eps_r = 1;
    auto ground = Medium(fd, eps_r, 1);
    auto lm = HalfSpace(fd, ground);
    EmMode mode = EmMode::TM;
    bool direct_term = true;

    VectorR3 r_ = {0, 0, 1};
    r_ *= fd.lambda_0;

    arma::vec x_vals = arma::linspace(-3, 3, 70) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-3, 3, 70) * fd.lambda_0;

    printf("x  z  g_re \n");
    for (auto &x : x_vals)
    {
        for (auto &z : z_vals)
        {
            VectorR3 r = {x, 0, z};
            auto g_lm = sgf::lm_generic_spatial(lm, r, r_, mode, direct_term);
            // auto g_0 = sgf::free_space(lm.media.back(), r, r_);

            printf("%.8f  %.8f  %.8f  \n", x, z, g_lm.real());
        }
        printf("\n");
    }

    return 0;
}
