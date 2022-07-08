#include <mthesis.hpp>

#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    FrequencyDomain fd(1e9);
    real w = 1.0 * fd.lambda_0;
    fmm::Params params(fd, w);

    // Source points: groups "far" above xy-plane.
    auto src_pts = fmm::rand_pts_in_group(params, {0, 0, 20}, 500);

    // Observation points: in xy-plane.
    arma::vec x_vals = arma::logspace(1, 4, 60) * fd.lambda_0;
    arma::vec y_vals = arma::logspace(1, 4, 60) * fd.lambda_0;
    std::vector<VectorR3> obs_pts;
    obs_pts.reserve(x_vals.n_elem * y_vals.n_elem);
    for (const auto &x : x_vals) {
        for (const auto &y : y_vals) {
            obs_pts.emplace_back(VectorR3 {x, y, 0.0});
        }
    }

    // Excitation vector.
    std::vector<cmplx> I(src_pts.size(), 1.0);

    // FMM solution.
    fmm::freespace::FMM my_fmm(params, src_pts, obs_pts);
    auto V_fmm = my_fmm.calc_product(I);

    // MoM solution.
    mom::MoM my_mom(fd, src_pts, obs_pts);
    auto V_mom = my_mom.calc_product(I);

    // Print informations into header comments.
    printf("# src groups:          %lu\n", my_fmm.src_groups.size());
    printf("# obs groups:          %lu\n", my_fmm.obs_groups.size());
    printf("# TOPs:                %lu\n", my_fmm.top_all.size());
    printf("# plane waves:         %lu\n", my_fmm.es.k_hat.size());
    printf("# multipole order: L = %u\n", my_fmm.L);
    // Print error.
    printf("x_by_lambda_0 y_by_lambda_0 rel_err_db\n");
    unsigned n = 0;
    for (const auto &x : x_vals) {
        for (const auto &y : y_vals) {
            real rel_err_db = calc_rel_err_db(V_fmm[n], V_mom[n]);
            n++;
            printf("%.6e %.6e %.6f\n",
                   x / fd.lambda_0,
                   y / fd.lambda_0,
                   rel_err_db);
        }
        printf("\n");
    }

    return 0;
}
