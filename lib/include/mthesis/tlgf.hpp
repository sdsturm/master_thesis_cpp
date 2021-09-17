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
               bool direct_term = true);

    const LayeredMedium &lm;
    const EmMode mode;
    const bool direct_term;
};

// Equivalent transmission line network Green's functions.
cmplx V_i(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet);
cmplx V_i(cmplx k_rho, real z, real z_, TLGFParams p);

cmplx I_i(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet);
cmplx I_i(cmplx k_rho, real z, real z_, TLGFParams p);

cmplx I_v(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet);
cmplx I_v(cmplx k_rho, real z, real z_, TLGFParams p);

cmplx V_v(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet);
cmplx V_v(cmplx k_rho, real z, real z_, TLGFParams p);

// Generalzied Sommerfeld identity.
cmplx generic_sgf(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet);
cmplx generic_sgf(cmplx k_rho, real z, real z_, TLGFParams p);

// *****************************************************************************

namespace utils {

// Basic quantities like k_z, impedance, reflection coefficients, ...
void calc_k_z_select_sheet(std::vector<cmplx> &k_z, RiemannSheet sheet);

std::vector<cmplx> calc_k_z(const TLGFParams &p,
                            cmplx k_rho,
                            RiemannSheet sheet);

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
    Internals(const TLGFParams &p,
              cmplx k_rho,
              bool dual_solution,
              RiemannSheet sheet);

    std::vector<cmplx> k_z;
    std::vector<cmplx> Z;
    std::vector<cmplx> theta;
    std::vector<cmplx> Gamma_u;
    std::vector<cmplx> Gamma_d;
};

struct LayerCoords
{
    LayerCoords(const LayeredMedium &lm, real z, real z_);

    real z;
    real z_;
    int m;
    int n;
};

// TLGF for z inside source layer.
std::vector<cmplx> calc_R(const Internals &d, const LayerCoords &c);

std::vector<cmplx> calc_zeta(const TLGFParams &p, const LayerCoords &c, real z);

cmplx calc_D(const TLGFParams &p, const LayerCoords &c, const Internals &d);

cmplx V_i_src_wo_prefac(const TLGFParams &p,
                        const LayerCoords &c,
                        real z,
                        const Internals &d);

cmplx V_i_src(const TLGFParams &p,
              const LayerCoords &c,
              real z,
              const Internals &d);

cmplx I_i_src(const TLGFParams &p,
              const LayerCoords &c,
              real z,
              const Internals &d);

// Transmission through the structure if n != m.
cmplx calc_tau_ud(int n,
                  const std::vector<cmplx> &Gamma_ud,
                  const std::vector<cmplx> &theta);

cmplx calc_tau_prod_d(const LayerCoords &c, const Internals &d);

cmplx calc_tau_prod_u(const LayerCoords &c, const Internals &d);

cmplx factor_eq_126(const TLGFParams &p,
                    const LayerCoords &c,
                    std::function<cmplx(cmplx, cmplx)> pm_operator,
                    const Internals &d);

cmplx factor_eq_117(const TLGFParams &p,
                    const LayerCoords &c,
                    std::function<cmplx(cmplx, cmplx)> pm_operator,
                    const Internals &d);

cmplx T_d(const TLGFParams &p,
          const LayerCoords &c,
          std::function<cmplx(cmplx, cmplx)> pm_operator,
          const Internals &d);

cmplx T_u(const TLGFParams &p,
          const LayerCoords &c,
          std::function<cmplx(cmplx, cmplx)> pm_operator,
          const Internals &d);

// General TLGF.
cmplx V_i_base_wo_prefac(const TLGFParams &p,
                         const LayerCoords &c,
                         const Internals &d);

cmplx V_i_base(const TLGFParams &p,
               const LayerCoords &c,
               const Internals &d);

cmplx I_i_base(const TLGFParams &p,
               const LayerCoords &c,
               const Internals &d);

} // namespace utils

} // namespace mthesis::tlgf

#endif
