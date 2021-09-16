#ifndef MTHESIS_TLGF_HPP
#define MTHESIS_TLGF_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

#include <functional>

// Reference: Michalski2005, Section 7.

namespace mthesis::tlgf {

// Equivalent transmission line network Green's functions.
cmplx V_i(const LayeredMedium &lm,
          cmplx k_rho,
          real z,
          real z_,
          EmMode type,
          bool direct_term = true,
          RiemannSheet sheet = RiemannSheet::I);

cmplx I_i(const LayeredMedium &lm,
          cmplx k_rho,
          real z,
          real z_,
          EmMode type,
          bool direct_term = true,
          RiemannSheet sheet = RiemannSheet::I);

cmplx I_v(const LayeredMedium &lm,
          cmplx k_rho,
          real z,
          real z_,
          EmMode type,
          bool direct_term = true,
          RiemannSheet sheet = RiemannSheet::I);

cmplx V_v(const LayeredMedium &lm,
          cmplx k_rho,
          real z,
          real z_,
          EmMode type,
          bool direct_term = true,
          RiemannSheet sheet = RiemannSheet::I);

// Hide implementation details in nested namespace utils.
namespace utils {

// Basic quantities like k_z, impedance, reflection coefficients, ...
void calc_k_z_select_sheet(std::vector<cmplx> &k_z, RiemannSheet sheet);

std::vector<cmplx> calc_k_z(const LayeredMedium &lm,
                            cmplx k_rho,
                            RiemannSheet sheet);

std::vector<cmplx> calc_Z(const LayeredMedium &lm,
                          const std::vector<cmplx> &k_z,
                          EmMode type);

std::vector<cmplx> calc_Y(const LayeredMedium &lm,
                          const std::vector<cmplx> &k_z,
                          EmMode type);

std::vector<cmplx> calc_theta(const LayeredMedium &lm,
                              const std::vector<cmplx> &k_z);

cmplx calc_Gamma_fresnel(const std::vector<cmplx> &Z, int i, int j);

std::vector<cmplx> calc_Gamma_u(const LayeredMedium &lm,
                                const std::vector<cmplx> &Z,
                                const std::vector<cmplx> &theta);

std::vector<cmplx> calc_Gamma_d(const LayeredMedium &lm,
                                const std::vector<cmplx> &Z,
                                const std::vector<cmplx> &theta);

struct Internals
{
    std::vector<cmplx> k_z;
    std::vector<cmplx> Z;
    std::vector<cmplx> theta;
    std::vector<cmplx> Gamma_u;
    std::vector<cmplx> Gamma_d;

    Internals(const LayeredMedium &lm,
              cmplx k_rho,
              EmMode mode,
              bool dual_solution,
              RiemannSheet riemann_sheet);
};

// TLGF for z inside source layer.
std::vector<cmplx> calc_R(const Internals &d, int n);

std::vector<cmplx> calc_zeta(const LayeredMedium &lm, real z, real z_, int n);

cmplx calc_D(const LayeredMedium &lm, int n, const Internals &d);

cmplx V_i_src(const LayeredMedium &lm,
              real z,
              real z_,
              int n,
              const Internals &d,
              bool direct_term);

cmplx I_i_src(const LayeredMedium &lm,
              real z,
              real z_,
              int n,
              const Internals &d,
              bool direct_term);

// Transmission through the structure if n != m.
cmplx calc_tau_ud(int n,
                  const std::vector<cmplx> &Gamma_ud,
                  const std::vector<cmplx> &theta);

cmplx calc_tau_prod_d(int m, int n, const Internals &d);

cmplx calc_tau_prod_u(int m, int n, const Internals &d);

cmplx factor_eq_126(const LayeredMedium &lm,
                    real z,
                    int m,
                    std::function<cmplx(cmplx, cmplx)> pm_operator,
                    const Internals &d);

cmplx factor_eq_117(const LayeredMedium &lm,
                    real z,
                    int m,
                    std::function<cmplx(cmplx, cmplx)> pm_operator,
                    const Internals &d);

cmplx T_d(const LayeredMedium &lm,
          real z,
          int m,
          int n,
          std::function<cmplx(cmplx, cmplx)> pm_operator,
          const Internals &d);

cmplx T_u(const LayeredMedium &lm,
          real z,
          int m,
          int n,
          std::function<cmplx(cmplx, cmplx)> pm_operator,
          const Internals &d);

// General TLGF.
cmplx V_i_base_generic(const LayeredMedium &lm,
                       real z,
                       real z_,
                       const Internals &d,
                       bool direct_term);

cmplx V_i_base(const LayeredMedium &lm,
               real z,
               real z_,
               const Internals &d,
               bool direct_term);

cmplx I_i_base(const LayeredMedium &lm,
               real z,
               real z_,
               const Internals &d,
               bool direct_term);

} // namespace utils

} // namespace mthesis::tlgf

#endif
