#ifndef MTHESIS_SI_NONSPECTRAL_HPP
#define MTHESIS_SI_NONSPECTRAL_HPP

#include <mthesis/si/sommerfeld_integral.hpp>

namespace mthesis::si::nonspectral {

// Hide implementation details in nested namespace.
namespace utils {

cmplx normalized_hankel(real n, cmplx z);

cmplx calc_folded_sgf(const SommerfeldIntegral &si, const LMCoords &c);

cmplx calc_F_eq81(cmplx s, const SommerfeldIntegral &si, const LMCoords &c);




// Hide implementation details in nested namespace.
namespace tm {

cmplx get_k_p(const LayeredMedium &lm);

cmplx get_s_p(const LayeredMedium &lm, cmplx k_p);

cmplx calc_k_p_residue(const SommerfeldIntegral &si, real z, real z_, cmplx k_p);

cmplx calc_numerical_distance(const LayeredMedium &lm, cmplx k_p, real rho);

cmplx calc_F_C8(cmplx p);

cmplx calc_B_p(cmplx k_p, cmplx s_p, cmplx R_p, real n, real rho);

cmplx integrand_I_p(cmplx s, cmplx s_p, cmplx B_p, real rho);

} // namespace tm




// Hide implementation details in nested namespace.
namespace te {

} // namespace te




} // namespace utils

} // namespace mthesis::si::nonspectral

#endif // MTHESIS_SI_NONSPECTRAL_HPP
