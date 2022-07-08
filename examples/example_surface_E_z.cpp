#include <mthesis.hpp>
#include "./misc/surface_E_z.hpp"

#include <armadillo>

#include "./../submodules/gnuplot-iostream/gnuplot-iostream.h"

using namespace mthesis;

int main()
{
    // Dry ground case from p.22 in Michalski2016b.
    FrequencyDomain fd(10e6);
    cmplx eps_r = cmplx_permittivity(fd, 3.0, 0.1e-3);
    Medium ground(fd, eps_r);
    HalfSpace lm(fd, ground);

    VectorR3 r_ = {0, 0, 1}; // In meter.

    arma::vec rho_vals = arma::logspace(-1, 3, 100) * fd.lambda_0;

    real h = 1.0;
    auto si = get_E_z_surf_VED_si(lm, h);

    arma::vec abs_E_z_VED(rho_vals.n_elem);
    using mthesis::si::axial_transmission::eval_si_along_sip;
    for (size_t n = 0; n < rho_vals.size(); n++) {
        VectorR3 r = {rho_vals(n), 0, 0};
        abs_E_z_VED[n] = abs(eval_si_along_sip(si, r, r_));
    }

    Gnuplot gp;
    gp << "set grid\n";
    gp << "set logscale xy\n";
    gp << "plot '-' with lines title '|E_z|'\n";
    gp.send1d(boost::make_tuple(rho_vals, abs_E_z_VED));

    return 0;
}
