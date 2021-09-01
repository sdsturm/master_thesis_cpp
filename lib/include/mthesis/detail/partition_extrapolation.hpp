#ifndef MTHESIS_PARTITION_EXTRAPOLATION_HPP
#define MTHESIS_PARTITION_EXTRAPOLATION_HPP

#include <mthesis/detail/definitions.hpp>
#include <mthesis/detail/spectral_gf.hpp>
#include <mthesis/detail/pe_params.hpp>

namespace mthesis
{
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
                          PEParams params);

    cmplx pe_levin_sidi(const SpectralGF &gf,
                        real nu,
                        real rho,
                        real a,
                        PEParams params);

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
                               PEParams params);

    cmplx pe_mosig_michalski(const SpectralGF &gf,
                             real nu,
                             real rho,
                             real a,
                             PEParams params);

} // namespace mthesis

#endif
