#include <mthesis.hpp>

#include <gsl/gsl_const_mksa.h>
#include <boost/timer/timer.hpp>

#include <cstdio>

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


    // Excitation vector.
    std::vector<cmplx> I(src_pts.size(), 1.0);

    boost::timer::auto_cpu_timer timer;

    // FMM solution.
    timer.start();
    fmm::FreeSpaceFMM my_fmm(params, src_pts, obs_pts);
    printf("FMM setup time: %s seconds\n", timer.format(9, "%u").c_str());

    timer.start();
    auto V_fmm = my_fmm.calc_product(I);
    printf("FMM solution time: %s seconds\n\n", timer.format(9, "%u").c_str());

    // MoM solution.
    timer.start();
    mom::MoM my_mom(fd, src_pts, obs_pts);
    printf("MoM setup time: %s seconds\n", timer.format(9, "%u").c_str());

    timer.start();
    auto V_mom = my_mom.calc_product(I);
    printf("MoM solution time: %s seconds\n\n", timer.format(9, "%u").c_str());

    // Compute error.
    std::vector<real> rel_err_db(obs_pts.size());
    for (size_t n = 0; n < obs_pts.size(); n++) {
        rel_err_db[n] = calc_rel_err_db(V_fmm[n], V_mom[n]);
    }

    printf("Maximum error: %.2f dB\n",
           *std::max_element(rel_err_db.begin(), rel_err_db.end()));

    return 0;
}
