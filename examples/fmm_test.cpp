#include <mthesis.hpp>

#include <gsl/gsl_const_mksa.h>

#include <iostream>

using namespace mthesis;

int main()
{
    FrequencyDomain fd(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    real w = fd.lambda_0;

    fmm::Params params(fd, w);

    unsigned n_pts = 10;

    auto src_pts = fmm::rand_pts_in_group(params, {0, 0, 4}, n_pts);
    fmm::append_pts(src_pts, fmm::rand_pts_in_group(params, {5, 0, 0}, n_pts));
    fmm::append_pts(src_pts, fmm::rand_pts_in_group(params, {-6, 3, -2}, n_pts));

    auto obs_pts = fmm::rand_pts_in_group(params, {7, 0, 10}, n_pts);
    fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, {3, 20, 15}, n_pts));
    fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, {30, 10, 5}, n_pts));
    fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, {-50, -10, -10}, n_pts));

    fmm::FreeSpaceFMM my_fmm(params, src_pts, obs_pts);

    mom::MoM my_mom(fd, src_pts, obs_pts);

    std::vector<cmplx> I(src_pts.size(), 1.0);

    auto V_fmm = my_fmm.calc_product(I);
    auto V_mom = my_mom.calc_product(I);

    for (size_t n = 0; n < obs_pts.size(); n++)
        std::cout << V_fmm[n] << "\t\t" << V_mom[n] << "\n";

    return 0;
}
