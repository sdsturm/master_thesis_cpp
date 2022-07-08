#ifndef MTHESIS_PE_PARTITION_EXTRAPOLATION_HPP
#define MTHESIS_PE_PARTITION_EXTRAPOLATION_HPP

#include <mthesis/definitions.hpp>

namespace mthesis::si::pe {

struct Params
{
    Params();
    void check() const;

    unsigned max_intervals;
    real tol;
};

real get_first_zero(real nu, real a, real rho);

cmplx levin_sidi(std::function<cmplx(real)> f,
                 real rho,
                 real a,
                 Params params);

cmplx mosig_michalski(std::function<cmplx(real)> f,
                      real alpha,
                      real zeta,
                      real rho,
                      real a,
                      Params params);

// Hide helper functions in nested namespace utils.
namespace utils {

real mc_mahon(real nu, unsigned m);

std::vector<real> get_j_table(real nu);

std::vector<real> get_xi(real a, real rho, unsigned max_intervals);

bool check_converged(cmplx val, const std::vector<cmplx> &old, real tol);

cmplx levin_sidi_extrap(int k,
                        cmplx s_k,
                        cmplx omega_k,
                        const std::vector<real> &xi,
                        std::vector<cmplx> &A,
                        std::vector<cmplx> &B);

cmplx mosig_michalski_extrap(real mu,
                             int k,
                             cmplx s_k,
                             real Omega_k,
                             const std::vector<real> &xi,
                             std::vector<cmplx> &R);

} // utils

} // namespace mthesis::si::pe

#endif
