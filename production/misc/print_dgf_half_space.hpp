#ifndef MTHESIS_PRINT_DGF_HALF_SPACE_HPP
#define MTHESIS_PRINT_DGF_HALF_SPACE_HPP

#include <mthesis.hpp>

#include <armadillo>
#include <gsl/gsl_const_mksa.h>

#include <cstdio>

namespace mthesis {

inline void print_dgf_half_space(const VectorR3 &J)
{
    using std::complex_literals::operator""i;

    // Dry ground example from Michalski2015.
    FrequencyDomain fd(10e6);
    auto eps_r = cmplx_permittivity(fd, 3.0, 0.1e-3);
    Medium ground(fd, eps_r);
    HalfSpace hs(fd, ground);

    VectorR3 r_ = {0, 0, fd.lambda_0};

    arma::vec x_vals = arma::linspace(0, 10, 70) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-5, 5, 70) * fd.lambda_0;

    printf("x/lambda_0 z/lambda_0, |Re(E)| Re(E_x) Re(E_y) Re(E_z)\n");
    for (const auto &x : x_vals) {
        for (const auto &z : z_vals) {
            VectorR3 r = {x, 0, z};
            DyadC3 G = gf::dyadic::layered_media::G_EJ(hs, r, r_);
            VectorC3 E = G * J;
            printf("%.6f %.6f %.6e %.6e %.6e %.6e\n",
                   x / fd.lambda_0,
                   z / fd.lambda_0,
                   arma::norm(arma::real(E)),
                   E(0).real(),
                   E(1).real(),
                   E(2).real());
        }
        printf("\n");
    }
}

} // namespace mthesis

#endif // MTHESIS_PRINT_DGF_HALF_SPACE_HPP
