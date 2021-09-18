#include <mthesis.hpp>

#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    using std::complex_literals::operator""i;

    // Dry ground example from Michalski2015.
    FrequencyDomain fd(10e6);
    auto eps_r = cmplx_permittivity(fd, 3.0, 0.1e-3);
    Medium ground(fd, eps_r);
    HalfSpace hs(fd, ground);

    VectorR3 r_ = {0, 0, 1.0 * fd.lambda_0};

    VectorR3 J_VED = {0, 0, 1};
    VectorR3 J_HED = {1, 0, 0};

    unsigned n_pts = 70;
    arma::vec x_vals = arma::linspace(-5, 5, n_pts) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-3, 7, n_pts) * fd.lambda_0;

    printf("x_by_lambda_0 z_by_lambda_0 abs_re_E_VED abs_re_E_HED\n");
    for (const auto &x : x_vals) {
        for (const auto &z : z_vals) {
            VectorR3 r = {x, 0, z};

            DyadC3 G_EJ = gf::dyadic::layered_media::G_EJ(hs, r, r_);
            VectorC3 E_VED = G_EJ * J_VED;
            VectorC3 E_HED = G_EJ * J_HED;

            printf("%.6e %.6e %.6e %.6e\n",
                   x / fd.lambda_0,
                   z / fd.lambda_0,
                   arma::norm(arma::real(E_VED)),
                   arma::norm(arma::real(E_HED)));
        }
        printf("\n");
    }

    return 0;
}
