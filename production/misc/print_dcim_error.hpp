#ifndef MTHESIS_PRINT_DCIM_ERROR_HPP
#define MTHESIS_PRINT_DCIM_ERROR_HPP

#include <mthesis.hpp>

#include <armadillo>

namespace mthesis {

inline void print_dcim_error(const arma::vec &rho_vals_by_lambda_0,
                             const arma::vec &z_vals_by_lambda_0)
{
    using std::complex_literals::operator""i;

    FrequencyDomain fd(1e9);
    Medium ground(fd, 3.0 - 0.2i, 1.0);
    HalfSpace hs(fd, ground);

    real z_ = 0.1 * fd.lambda_0;
    VectorR3 r_ = {0, 0, z_};

    arma::vec rho_vals = rho_vals_by_lambda_0 * fd.lambda_0;
    arma::vec z_vals = z_vals_by_lambda_0 * fd.lambda_0;

    bool direct_term = false;
    EmMode mode = EmMode::TM; // Sommerfeld pole exists.
    real nu = 0;

    auto si = gf::scalar::layered_media::get_sommerfeld_integral(
                hs, nu, mode, direct_term);

    auto dcim_3lv2 = dcim::ThreeLevelV2(si);

    real rel_err_db_max = -INFINITY;
    real rho_max_err, z_max_err;
    printf("rho/lambda_0 z/lambda_0 rel_err_db\n");
    for (const auto &z : z_vals) {
        auto ce_vecs_levels = dcim_3lv2.get_exponentials(z, z_);
        for (const auto &rho : rho_vals) {
            auto val_ref = si.eval_si_along_sip(rho, z, z_);
            auto val_dcim_3lv2 = dcim_3lv2.get_spatial_gf(ce_vecs_levels, rho);

            auto rel_err_db = calc_rel_err_db(val_dcim_3lv2, val_ref);

            if (rel_err_db > rel_err_db_max) {
                rel_err_db_max = rel_err_db;
                rho_max_err = rho;
                z_max_err = z;
            }

            printf("%.6f %.6f %.6f\n",
                   rho / fd.lambda_0,
                   z / fd.lambda_0,
                   rel_err_db);
        }
        printf("\n");
    }

//    printf("Maximum error %.2f dB at rho = %.4e*lambda_0, z = %.4e*lambda_0\n",
//           rel_err_db_max, rho_max_err / fd.lambda_0, z_max_err / fd.lambda_0);
}

}

#endif // MTHESIS_PRINT_DCIM_ERROR_HPP
