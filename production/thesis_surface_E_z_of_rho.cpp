#include <mthesis.hpp>

#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    // Dry ground case from p.22 in Michalski2016b.
    FrequencyDomain fd(10e6);
    cmplx eps_r = cmplx_permittivity(fd, 3.0, 0.1e-3);
    Medium ground(fd, eps_r);
    HalfSpace lm(fd, ground);

    VectorR3 r_ = {0, 0, 1}; // In meter.

    arma::vec rho_vals = arma::logspace(-1, 3, 150) * fd.lambda_0;

    printf("rho_by_lambda_0 abs_E_z\n");
    for (const auto &rho : rho_vals) {
        VectorR3 r = {rho, 0, 0};
        auto G_EJ = gf::dyadic::layered_media::G_EJ(lm, r, r_);
        auto E_z_VED = G_EJ(2, 2);

        printf("%.6e %.6e\n", rho / fd.lambda_0, abs(E_z_VED));
    }
    return 0;
}
