#include <mthesis.hpp>

#include <gsl/gsl_const_mksa.h>
#include <armadillo>

#include <filesystem>
#include <iostream>

#include "./../submodules/gnuplot-iostream/gnuplot-iostream.h"

using namespace mthesis;

int main()
{
    using std::complex_literals::operator""i;

    // Problem specification.
    auto fd = FrequencyDomain(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    Medium medium(fd, 3.0 - 0.2i, 1);
    std::vector<real> z_interfces = {-INFINITY,
                                     -fd.lambda_0,
                                     fd.lambda_0,
                                     INFINITY};
    std::vector<Medium> media = {Vacuum(fd), medium, Vacuum(fd)};
    LayeredMedium lm(fd, z_interfces, media, BC::open, BC::open);

    real k_rho = 0.0;
    EmMode type = EmMode::TM;

    real z_ = 2.0 * fd.lambda_0;

    arma::vec z = arma::linspace(-3, 3, 1e3) * fd.lambda_0;

    // TLGF computation.
    arma::cx_vec V(z.n_elem), I(z.n_elem), Z_by_Z_0(z.n_elem);

    const real Z_0 = sqrt(GSL_CONST_MKSA_VACUUM_PERMEABILITY /
                          GSL_CONST_MKSA_VACUUM_PERMITTIVITY);

    using namespace mthesis::tlgf;
    for (arma::uword n = 0; n < z.n_elem; n++)
    {
        V(n) = V_i(k_rho, z(n), z_, TLGFParams(lm, type));
        I(n) = I_i(k_rho, z(n), z_, TLGFParams(lm, type));
        Z_by_Z_0(n) = V(n) / I(n) / Z_0;
    }

    // Plot results.
    Gnuplot gp;
    gp << "set multiplot layout 3,1\n";
    gp << "set grid\n";
    gp << "plot '-' with lines title 'Re(V)', '-' with lines title 'Im(V)'\n";
    gp.send1d(boost::make_tuple(z, arma::vec(arma::real(V))));
    gp.send1d(boost::make_tuple(z, arma::vec(arma::imag(V))));
    gp << "plot '-' with lines title 'Re(I)', '-' with lines title 'Im(I)'\n";
    gp.send1d(boost::make_tuple(z, arma::vec(arma::real(I))));
    gp.send1d(boost::make_tuple(z, arma::vec(arma::imag(I))));
    gp << "plot '-' with lines title 'Re(Z/Z_0)', '-' with lines title 'Im(Z/Z_0)'\n";
    gp.send1d(boost::make_tuple(z, arma::vec(arma::real(Z_by_Z_0))));
    gp.send1d(boost::make_tuple(z, arma::vec(arma::imag(Z_by_Z_0))));

    return 0;
}
