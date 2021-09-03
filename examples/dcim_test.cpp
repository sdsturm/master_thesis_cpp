#include <mthesis.hpp>

#include <gsl/gsl_const_mksa.h>

using namespace mthesis;

VectorR3 rand_point(real factor)
{
    return (2 * arma::randu(3) - 1) * factor;
}

int main()
{
    using std::complex_literals::operator""i;

    auto fd = FrequencyDomain(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    auto ground = Medium(fd, 3.0, 1.0);
    auto lm = HalfSpace(fd, ground);
    auto mode = EmMode::TM;

     auto r = rand_point(5 * fd.lambda_0);
     auto r_ = rand_point(5 * fd.lambda_0);
//    VectorR3 r = {5, 0, 1};
//    VectorR3 r_ = {0, 0, 1};
    r *= fd.lambda_0;
    r_ *= fd.lambda_0;

    std::cout << "r = \n" << r << "\n";
    std::cout << "r_ = \n" << r_ << "\n";

    bool direct_term = false;
    real nu = 0;
    
    auto gf = sgf::lm_get_generic_spectral_gf(lm, r, r_, mode, direct_term);

    auto coords = LayeredMediumCoords(r, r_);
    auto val_ref = si::eval_along_sip(gf, nu, coords.rho);

    auto ce_levels = dcim::threelevelv2::three_level_v2(gf);

    auto val_dcim = dcim::utils::get_spectral_gf(ce_levels, coords.rho, fd.k_0);

    std::cout << "val_ref:        " << val_ref << "\n";
    std::cout << "val_dcim:       " << val_dcim << "\n";
    std::cout << "Relative error: " << calc_rel_err_db(val_dcim, val_ref) << " dB\n";

    return 0;
}
