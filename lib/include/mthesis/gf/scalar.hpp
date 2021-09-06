#ifndef MTHESIS_GF_SCALAR_HPP
#define MTHESIS_GF_SCALAR_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/sommerfeld_integrals.hpp>

namespace mthesis::gf::scalar {

namespace free_space {

cmplx G_0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

cmplx G_0(const Medium &medium, const VectorR3 &r, const VectorC3 &r_);

} // namespace free_space

namespace layered_media {

cmplx generic_spectral(const LayeredMedium &lm,
                       real z,
                       real z_,
                       cmplx k_rho,
                       EmMode mode,
                       bool direct_term);

SommerfeldIntegral get_sommerfeld_integral(const LayeredMedium &lm,
                                           real nu,
                                           EmMode mode,
                                           bool direct_term,
                                           SiParams si_params = SiParams());

} // namespace layered_media

} // namespace mthesis::gf::scalar

#endif // MTHESIS_GF_SCALAR_HPP
