#include <mthesis.hpp>

#include <gsl/gsl_const_mksa.h>

using namespace mthesis;

int main()
{
    FrequencyDomain fd(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    real w = fd.lambda_0;

    fmm::Params params(fd, w);

    unsigned n_pts = 300;

    auto src_pts = fmm::rand_pts_in_group(params, {0, 0, 4}, n_pts);
    fmm::append_pts(src_pts, fmm::rand_pts_in_group(params, {5, 0, 0}, n_pts));
    fmm::append_pts(src_pts, fmm::rand_pts_in_group(params, {-6, 3, -2}, n_pts));

    auto obs_pts = fmm::rand_pts_in_group(params, {7, 0, 10}, n_pts);
    fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, {3, 20, 15}, n_pts));
    fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, {30, 10, 5}, n_pts));
    fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, {-50, -10, -10}, n_pts));

    fmm::FreeSpaceFMM my_fmm(params, src_pts, obs_pts);


    return 0;
}
