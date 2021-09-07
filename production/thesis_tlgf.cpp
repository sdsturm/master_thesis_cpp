#include <mthesis.hpp>

#include <armadillo>
#include <gsl/gsl_const_mksa.h>

#include <cstdio>

using namespace mthesis;

int main()
{
    using std::complex_literals::operator""i;

    FrequencyDomain fd(1e9);
    Medium ground(fd, 9.0 - 1.0i, 1);
    HalfSpace lm(fd, ground);

    real z_ = 1.0 * fd.lambda_0;

    arma::vec z_vals = arma::linspace(-2, 4, 800) * fd.lambda_0;

    auto type = EmMode::TM;
    bool direct_term = true;
    real k_rho = 0;
    real Z_0 = std::sqrt(GSL_CONST_MKSA_VACUUM_PERMEABILITY /
                         GSL_CONST_MKSA_VACUUM_PERMITTIVITY);

    printf("z_by_lambda_0 V_re V_im I_re I_im Z_rel_re Z_rel_im\n");
    for (const auto &z : z_vals) {
        cmplx V = tlgf::V_i(lm, k_rho, z, z_, type, direct_term);
        cmplx I = tlgf::I_i(lm, k_rho, z, z_, type, direct_term);
        cmplx Z_rel = V / I / Z_0;
        printf("%.6f %.6e %.6e %.6e %.6e %.6e %.6e\n",
               z / fd.lambda_0,
               V.real(), V.imag(),
               I.real(), I.imag(),
               Z_rel.real(), Z_rel.imag());


    }



    return 0;
}
