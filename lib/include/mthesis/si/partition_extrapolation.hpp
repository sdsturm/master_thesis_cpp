#ifndef MTHESIS_PE_PARTITION_EXTRAPOLATION_HPP
#define MTHESIS_PE_PARTITION_EXTRAPOLATION_HPP

#include <mthesis/definitions.hpp>

namespace mthesis::pe {

class Params
{
private:
    unsigned max_intervals;
    double tol;

public:
    Params() : max_intervals(15), tol(1e-8) {}

    void set_max_intervals(int max_intervals);
    void set_tol(double tol);

    unsigned get_max_intervals() const;
    double get_tol() const;
};


cmplx levin_sidi(std::function<cmplx(real)> f,
                 real nu,
                 real rho,
                 real a,
                 Params params);

cmplx mosig_michalski(std::function<cmplx(real)> f,
                      real alpha,
                      real zeta,
                      real nu,
                      real rho,
                      real a,
                      Params params);

// Hide helper functions in nested namespace utils.
namespace utils {

double mc_mahon(double nu, unsigned m);

double get_first_zero(double nu, double a, double rho);

cmplx integrate_gap(std::function<cmplx(real)> f,
                    double a,
                    double a_pe_start);

std::vector<double> get_xi(double a, double rho, unsigned max_intervals);

bool check_converged(cmplx val, const std::vector<cmplx> &old, double tol);

cmplx levin_sidi_extrap(int k,
                        cmplx s_k,
                        cmplx omega_k,
                        const std::vector<double> &xi,
                        std::vector<cmplx> &A,
                        std::vector<cmplx> &B);

cmplx levin_sidi_core(std::function<cmplx(real)> f,
                      double rho,
                      double a,
                      Params params);

cmplx mosig_michalski_extrap(double mu,
                             int k,
                             cmplx s_k,
                             double Omega_k,
                             const std::vector<double> &xi,
                             std::vector<cmplx> &R);

cmplx mosig_michalski_core(std::function<cmplx(real)> f,
                           double rho,
                           double a,
                           double alpha,
                           double zeta,
                           Params params);

} // namespace mthesis::pe::utils

} // namespace mthesis::pe

#endif