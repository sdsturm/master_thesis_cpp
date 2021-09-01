#include <mthesis/mthesis.hpp>

#include <gsl/gsl_const_mksa.h>

#include <iostream>

using namespace mthesis;

int main()
{
    // Problem specification.
    auto fd = FrequencyDomain(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    auto eps_r = cmplx_permittivity(fd, 3, 10e-3);
    auto ground = Medium(fd, eps_r, 1);
    auto lm = HalfSpace(fd, ground);
    EmMode mode = EmMode::TM;
    bool direct_term = true;

    VectorR3 r = arma::randu(3) * fd.lambda_0;
    VectorR3 r_ = arma::randu(3) * fd.lambda_0;

    auto g_lm = sgf::lm_generic_spatial(lm, r, r_, mode, direct_term);


    return 0;
}
