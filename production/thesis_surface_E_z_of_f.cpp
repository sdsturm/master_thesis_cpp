#include <mthesis.hpp>

#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    arma::vec f_vals = arma::logspace(5, 8, 600);

    // Distances in meter.
    VectorR3 r_ = {0, 0, 0};
    VectorR3 r = {50, 0, 0};

    printf("f abs_E_z\n");
    for (const auto &f : f_vals) {
        FrequencyDomain fd(f);
        cmplx eps_r = cmplx_permittivity(fd, 3, 0.1e-3);
        Medium ground(fd, eps_r);
        HalfSpace lm(fd, ground);
        auto G_EJ = gf::dyadic::layered_media::G_EJ(lm, r, r_);
        auto E_z_VED = G_EJ(2, 2);

        printf("%.6e %.6e\n", f, abs(E_z_VED));
    }

    return 0;
}
