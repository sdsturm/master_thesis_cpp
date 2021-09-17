#ifndef MTHESIS_TLGF_HPP
#define MTHESIS_TLGF_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

#include <functional>

// Reference: Michalski2005, Section 7.

namespace mthesis::tlgf {

struct TLGFParams
{
    TLGFParams(const LayeredMedium &lm,
               EmMode mode,
               bool direct_term = true,
               RiemannSheet sheet = RiemannSheet::I);

    const LayeredMedium &lm;
    const EmMode mode;
    const bool direct_term;
    RiemannSheet sheet;
};

// Equivalent transmission line network Green's functions.
cmplx V_i(cmplx k_rho, real z, real z_, TLGFParams p);

cmplx I_i(cmplx k_rho, real z, real z_, TLGFParams p);

cmplx I_v(cmplx k_rho, real z, real z_, TLGFParams p);

cmplx V_v(cmplx k_rho, real z, real z_, TLGFParams p);

cmplx generic_sgf(cmplx k_rho, real z, real z_, TLGFParams p);

// Hide implementation details in nested namespace utils.
namespace utils {

// Basic quantities like k_z, impedance, reflection coefficients, ...
void calc_k_z_select_sheet(std::vector<cmplx> &k_z, const TLGFParams &p);

std::vector<cmplx> calc_k_z(const TLGFParams &p, cmplx k_rho);

std::vector<cmplx> calc_Z(const TLGFParams &p, const std::vector<cmplx> &k_z);

std::vector<cmplx> calc_Y(const TLGFParams &p, const std::vector<cmplx> &k_z);

std::vector<cmplx> calc_theta(const TLGFParams &p,
                              const std::vector<cmplx> &k_z);

cmplx calc_Gamma_fresnel(const std::vector<cmplx> &Z, int i, int j);

std::vector<cmplx> calc_Gamma_u(const TLGFParams &p,
                                const std::vector<cmplx> &Z,
                                const std::vector<cmplx> &theta);

std::vector<cmplx> calc_Gamma_d(const TLGFParams &p,
                                const std::vector<cmplx> &Z,
                                const std::vector<cmplx> &theta);

struct Internals
{
    std::vector<cmplx> k_z;
    std::vector<cmplx> Z;
    std::vector<cmplx> theta;
    std::vector<cmplx> Gamma_u;
    std::vector<cmplx> Gamma_d;

    Internals(const TLGFParams &p, cmplx k_rho, bool dual_solution);
};

// TLGF for z inside source layer.
std::vector<cmplx> calc_R(const Internals &d, int n);

std::vector<cmplx> calc_zeta(const LayeredMedium &lm, real z, real z_, int n);

cmplx calc_D(const LayeredMedium &lm, int n, const Internals &d);

cmplx V_i_src_wo_prefac(const LayeredMedium &lm,
                        real z,
                        real z_,
                        int n,
                        const Internals &d,
                        bool direct_term);

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
cmplx V_i_base_wo_prefac(const LayeredMedium &lm,
                         real z,
                         real z_,
                         int m,
                         int n,
                         const Internals &d,
                         bool direct_term);

cmplx V_i_base(const LayeredMedium &lm,
               real z,
               real z_,
               int m,
               int n,
               const Internals &d,
               bool direct_term);

cmplx I_i_base(const LayeredMedium &lm,
               real z,
               real z_,
               int m,
               int n,
               const Internals &d,
               bool direct_term);

} // namespace utils

} // namespace mthesis::tlgf

#endif
