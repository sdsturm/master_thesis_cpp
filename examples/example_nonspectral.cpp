#include <mthesis.hpp>

#include <armadillo>

#include <iostream>

using namespace mthesis;

int main()
{
    FrequencyDomain fd(10e6);
    cmplx eps_r = cmplx_permittivity(fd, 3.0, 0.1e-3);
    Medium ground(fd, eps_r);
    HalfSpace lm(fd, ground);
    real nu = 0;
    EmMode mode = EmMode::TM;

    auto si = gf::scalar::layered_media::get_sommerfeld_integral(lm, nu, mode, true);

    arma::vec rho_vals = arma::logspace(-1, 3, 30) * fd.lambda_0;

    VectorR3 r_ = {0, 0, 1};	// In meter.

    std::cout << "val_ax val_ns rel_err_db\n";
    for (const auto &rho : rho_vals) {
        VectorR3 r = {rho, 0, 0};
        cmplx val_ax = si::axial_transmission::eval_si_along_sip(si, r, r_);


        LMCoords c(r, r_);
        cmplx val_ns = si::nonspectral::eval_nonspectral(si, c, mode);

        std::cout << val_ax << " " <<
                     val_ns << " " <<
                     calc_rel_err_db(val_ns, val_ax) <<
                     "\n";
    }

    return 0;
}
