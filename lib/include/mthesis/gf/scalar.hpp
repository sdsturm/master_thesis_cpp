#ifndef MTHESIS_GF_SCALAR_HPP
#define MTHESIS_GF_SCALAR_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/si/sommerfeld_integral.hpp>

namespace mthesis::gf::scalar {

namespace free_space {

cmplx G_0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

cmplx G_0(const Medium &medium, const VectorR3 &r, const VectorC3 &r_);

} // namespace free_space

namespace layered_media {

// Return a SommerfeldIntegral object which can be evaluated in spatial domain.
SommerfeldIntegral get_sommerfeld_integral(const LayeredMedium &lm,
                                           real nu,
                                           EmMode mode,
                                           bool direct_term);

} // namespace layered_media

} // namespace mthesis::gf::scalar

#endif // MTHESIS_GF_SCALAR_HPP
