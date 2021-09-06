#include <mthesis/fmm/gegenbauer.hpp>
#include <mthesis/gf.hpp>

#include <gsl/gsl_const_mksa.h>

#include <iostream>

using namespace mthesis;

int main()
{
    FrequencyDomain fd(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    Medium medium = Vacuum(fd);

    fmm::Params params(fd, 1.0 * fd.lambda_0);

    VectorR3 r =  fmm::rand_point_in_group(params, {10, 0, 0});
    VectorR3 r_ = fmm::rand_point_in_group(params, {0, 0, 0});

    cmplx ref = gf::scalar::free_space::G_0(medium, r, r_);

    unsigned L = 20;
    cmplx num = fmm::compute_gegenbauer_pw(params, r, r_, L);

    std::cout << "ref = " << ref << "\n";
    std::cout << "num = " << num << "\n";
    std::cout << "err = " << calc_rel_err_db(num, ref) << " dB\n";

    return 0;
}
