#ifndef MTHESIS_TLGF_HPP
#define MTHESIS_TLGF_HPP

#include <mthesis/detail/definitions.hpp>
#include <mthesis/detail/solution_domain.hpp>

#include <functional>

// Reference: Michalski2005, Section 7.

namespace mthesis
{
    enum class RiemannSheet
    {
        I,
        II,
        III,
        IV
    };

    void calc_k_z_select_sheet(std::vector<cmplx> &k_z,
                               RiemannSheet sheet);

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
                  RiemannSheet riemann_sheet = RiemannSheet::I);
    };

    std::vector<cmplx> calc_R(const Internals &d, int n);

    std::vector<cmplx> calc_zeta(const LayeredMedium &lm,
                                 real z,
                                 real z_,
                                 int n);

    cmplx calc_D(const LayeredMedium &lm,
                 int n,
                 const Internals &d);

    cmplx V_i_src_layer(const LayeredMedium &lm,
                        real z,
                        real z_,
                        int n,
                        const Internals &d);

    cmplx I_i_src_layer(const LayeredMedium &lm,
                        real z,
                        real z_,
                        int n,
                        const Internals &d);

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

    cmplx V_i_base(const LayeredMedium &lm,
                   real z,
                   real z_,
                   const Internals &d);

    cmplx I_i_base(const LayeredMedium &lm,
                   real z,
                   real z_,
                   const Internals &d);

    cmplx V_i(const LayeredMedium &lm,
              cmplx k_ρ,
              real z,
              real z_,
              EmMode type);

    cmplx I_i(const LayeredMedium &lm,
              cmplx k_ρ,
              real z,
              real z_,
              EmMode type);

    cmplx I_v(const LayeredMedium &lm,
              cmplx k_ρ,
              real z,
              real z_,
              EmMode type);

    cmplx V_v(const LayeredMedium &lm,
              cmplx k_ρ,
              real z,
              real z_,
              EmMode type);

} // namespace mthesis

#endif
