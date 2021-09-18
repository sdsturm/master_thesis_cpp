#include <mthesis.hpp>

#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    // Dry ground example from p. 22 in Michalski2016b.
    FrequencyDomain fd(10e6);
    cmplx eps_r = cmplx_permittivity(fd, 3, 0.1e-3);
    Medium ground(fd, eps_r);
    HalfSpace lm(fd, ground);
    real nu = 0;

    VectorR3 r_ = {0, 0, 1.0 * fd.lambda_0};

    arma::vec x_vals = arma::linspace(0, 10, 70) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-3, 7, 70) * fd.lambda_0;

    using mthesis::gf::scalar::layered_media::get_sommerfeld_integral;
    auto si_tm_full = get_sommerfeld_integral(lm, nu, EmMode::TM, true);
    auto si_tm_refl = get_sommerfeld_integral(lm, nu, EmMode::TM, false);
    auto si_te_full = get_sommerfeld_integral(lm, nu, EmMode::TE, true);
    auto si_te_refl = get_sommerfeld_integral(lm, nu, EmMode::TE, false);

    printf("x_by_lambda_0 z_by_lambda_0 g_tm_full_log g_tm_refl_re g_te_full_log g_te_refl_re\n");
    for (const auto &x : x_vals) {
        for (const auto &z : z_vals) {
            VectorR3 r = {x, 0, z};
            using mthesis::si::axial_transmission::eval_si_along_sip;
            cmplx g_tm_full = eval_si_along_sip(si_tm_full, r, r_);
            cmplx g_tm_refl = eval_si_along_sip(si_tm_refl, r, r_);
            cmplx g_te_full = eval_si_along_sip(si_te_full, r, r_);
            cmplx g_te_refl = eval_si_along_sip(si_te_refl, r, r_);

            printf("%.6e %.6e %.6e %.6e %.6e %.6e\n",
                   x / fd.lambda_0,
                   z / fd.lambda_0,
                   20 * log10(abs(g_tm_full)),
                   g_tm_refl.real(),
                   20 * log10(abs(g_te_full)),
                   g_te_refl.real());
        }
        printf("\n");
    }

    return 0;
}
