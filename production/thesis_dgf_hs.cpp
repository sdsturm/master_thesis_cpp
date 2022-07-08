#include <mthesis.hpp>

#include <armadillo>

#include <iostream>
#include <string>
#include <cassert>

std::string print_comp(unsigned i)
{
    switch (i) {
    case 0:
        return "x";
        break;
    case 1:
        return "y";
        break;
    case 2:
        return "z";
        break;
    default:
        assert(false);
    }
}

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

    unsigned n_pts = 60;
    arma::vec x_vals = arma::linspace(-5, 5, n_pts) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-3, 7, n_pts) * fd.lambda_0;

    // Print header.
    std::cout << "x_by_lambda_0 z_by_lambda_0 ";
    for (arma::uword r = 0; r < 3; r++) {
        for (arma::uword c = 0; c < 3; c++) {
            std::cout << "G_EJ_" << print_comp(r) << print_comp(c) << " ";
            std::cout << "G_HJ_" << print_comp(r) << print_comp(c) << " ";
        }
    }
    std::cout << "\n";

    // Run computations and print result line.
    for (const auto &x : x_vals) {
        for (const auto &z : z_vals) {
            VectorR3 r = {x, 0, z};

            DyadC3 G_EJ = gf::dyadic::layered_media::G_EJ(hs, r, r_);
            DyadC3 G_HJ = gf::dyadic::layered_media::G_HJ(hs, r, r_);

            std::cout << x / fd.lambda_0 << " " << z / fd.lambda_0 << " ";
            for (arma::uword r = 0; r < 3; r++) {
                for (arma::uword c = 0; c < 3; c++) {
                    std::cout << G_EJ(r, c).real() << " ";
                    std::cout << G_HJ(r, c).real() << " ";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    return 0;
}
