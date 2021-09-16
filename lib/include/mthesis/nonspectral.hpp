#ifndef MTHESIS_NONSPECTRAL_HPP
#define MTHESIS_NONSPECTRAL_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

namespace mthesis::nonspectral {

cmplx get_sommerfeld_pole(const LayeredMedium &lm, EmMode mode);

cmplx normalized_hankel_fun(real n, cmplx z);

cmplx faddeeva(cmplx z);

cmplx func_C8(cmplx p);

cmplx calc_I_q(const LayeredMedium &lm,
               cmplx k_p,
               real rho);

} // namespace mthesis::nonspectral

#endif // MTHESIS_NONSPECTRAL_HPP
