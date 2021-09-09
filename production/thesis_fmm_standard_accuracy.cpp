#include <mthesis.hpp>

#include <armadillo>
#include <gsl/gsl_const_mksa.h>
#include <boost/timer/timer.hpp>

#include <random>
#include <cstdio>

using namespace mthesis;

fmm::multiindex rand_mi()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 100);

    return fmm::multiindex {dist(gen), dist(gen), dist(gen)};
}

int main()
{
    FrequencyDomain fd(GSL_CONST_MKSA_SPEED_OF_LIGHT);
    real w = 1.0 * fd.lambda_0;
    fmm::Params params(fd, w);

    // Source points: groups "far" above xy-plane.
    auto src_pts = fmm::rand_pts_in_group(params, {0, 0, 20}, 500);

    // Observation points: in xy-plane.
    arma::vec x_vals = arma::linspace(-100, 100, 50) * w;
    arma::vec y_vals = arma::linspace(-100, 100, 50) * w;
    std::vector<VectorR3> obs_pts;
    obs_pts.reserve(x_vals.n_elem * y_vals.n_elem);
    for (const auto &x : x_vals) {
        for (const auto &y : y_vals) {
            obs_pts.emplace_back(VectorR3 {x, y, 0.0});
        }
    }

    // Excitation vector.
    std::vector<cmplx> I(src_pts.size(), 1.0);

    boost::timer::auto_cpu_timer timer;

    // FMM solution.
    printf("FMM setup: ");
    timer.start();
    fmm::freespace::FMM my_fmm(params, src_pts, obs_pts);
    printf("%s seconds\n", timer.format(9, "%u").c_str());

    printf("FMM solution: ");
    timer.start();
    auto V_fmm = my_fmm.calc_product(I);
    printf("%s seconds\n", timer.format(9, "%u").c_str());

    // MoM solution.
    printf("MoM setup: ");
    timer.start();
    mom::MoM my_mom(fd, src_pts, obs_pts);
    printf("%s\n", timer.format(9, "%u").c_str());

    printf("MoM solution: ");
    timer.start();
    auto V_mom = my_mom.calc_product(I);
    printf("%s\n", timer.format(9, "%u").c_str());

    // Compute error.
    std::vector<real> rel_err_db(obs_pts.size());
    for (size_t n = 0; n < obs_pts.size(); n++) {
        rel_err_db[n] = calc_rel_err_db(V_fmm[n], V_mom[n]);
    }

    printf("Maximum error: %.2f dB\n",
           *std::max_element(rel_err_db.begin(), rel_err_db.end()));

    return 0;
}
