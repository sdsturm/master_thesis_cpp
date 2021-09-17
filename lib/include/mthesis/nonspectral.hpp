#ifndef MTHESIS_NONSPECTRAL_HPP
#define MTHESIS_NONSPECTRAL_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/sommerfeld_integrals.hpp>

namespace mthesis::nonspectral {


cmplx normalized_hankel_fun(real n, cmplx z);

cmplx calc_f_folded(cmplx k_rho,
                    const std::function<cmplx (cmplx, RiemannSheet)> &f,
                    const LayeredMedium &lm);

cmplx calc_F_eq81(cmplx s,
                  int n,
                  const LayeredMedium &lm,
                  const std::function<cmplx (cmplx, RiemannSheet)> &f,
                  real rho);


cmplx eval_nonspectral(const LayeredMedium &lm,
                       const std::function<cmplx (cmplx, RiemannSheet)> &f,
                       real rho);

namespace tm_case {

cmplx get_sommerfeld_pole(const LayeredMedium &lm, EmMode mode);

cmplx calc_s_p(const LayeredMedium &lm, cmplx k_p);

cmplx integrate_k_p_residue(std::function<cmplx (cmplx)> f, cmplx k_p);

cmplx cmplx_calc_num_dist_eqC5(const LayeredMedium &lm, cmplx k_p, real rho);

cmplx calc_F_eqC8(cmplx p);

cmplx calc_I_q(const LayeredMedium &lm,
               cmplx k_p,
               real rho);

cmplx calc_B_p(const LayeredMedium &lm,
               real n,
               std::function<cmplx (cmplx)> f,
               cmplx k_p,
               real rho);

cmplx calc_I_p(int n,
               const LayeredMedium &lm,
               const std::function<cmplx (cmplx, RiemannSheet)> &f,
               real rho,
               cmplx s_p,
               cmplx B_p);

cmplx calc_I_1_n(const LayeredMedium &lm,
                 int n,
                 const std::function<cmplx (cmplx, RiemannSheet)> &f,
                 real rho);

} // namespace tm_case

namespace te_case {

} // namespace te_case

} // namespace mthesis::nonspectral

#endif // MTHESIS_NONSPECTRAL_HPP
