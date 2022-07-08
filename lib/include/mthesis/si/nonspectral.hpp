#ifndef MTHESIS_SI_NONSPECTRAL_HPP
#define MTHESIS_SI_NONSPECTRAL_HPP

#include <mthesis/si/sommerfeld_integral.hpp>

namespace mthesis::si::nonspectral {

cmplx eval_nonspectral(const SommerfeldIntegral &si,
                       const LMCoords &c,
                       EmMode mode);

// =============================================================================

namespace utils {

cmplx normalized_hankel(real n, cmplx z);

cmplx calc_folded_sgf(const SommerfeldIntegral &si,
                      const LMCoords &c,
                      cmplx k_rho);

#if 0
cmplx calc_F_eq81(cmplx s, const SommerfeldIntegral &si, const LMCoords &c);
#endif

// =============================================================================

namespace tm {

cmplx calc_k_p(const LayeredMedium &lm);

cmplx calc_s_p(const LayeredMedium &lm, cmplx k_p);

#if 0
cmplx calc_R_p_quad(const SommerfeldIntegral &si, const LMCoords &c, cmplx k_p);
#endif

cmplx calc_R_p_closed_form(const LayeredMedium &lm, const LMCoords &c);

cmplx calc_numerical_distance(const LayeredMedium &lm, cmplx k_p, real rho);

cmplx calc_F_C8(cmplx p);

cmplx calc_B_p(cmplx k_p, cmplx s_p, cmplx R_p, real n, real rho);

cmplx integrand_I_p(const SommerfeldIntegral &si, const LMCoords &c, cmplx s,
                    cmplx s_p, cmplx B_p);

cmplx calc_I_p(const SommerfeldIntegral &si, const LMCoords &c,
               cmplx B_p, cmplx s_p);

cmplx calc_I1_eq85(const SommerfeldIntegral &si, const LMCoords &c);

} // namespace tm

// =============================================================================

namespace te {

#if 0
cmplx integrand_I_p(const SommerfeldIntegral &si, const LMCoords &c, cmplx s);

cmplx calc_I_p(const SommerfeldIntegral &si, const LMCoords &c);

cmplx calc_I1_eq85(const SommerfeldIntegral &si, const LMCoords &c);
#endif

} // namespace te

} // namespace utils

} // namespace mthesis::si::nonspectral

#endif // MTHESIS_SI_NONSPECTRAL_HPP
