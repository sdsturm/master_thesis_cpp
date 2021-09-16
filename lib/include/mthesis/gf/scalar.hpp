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

// Generic spectral Green's functions (generalized Sommerfeld identity).
cmplx generic_spectral_gf(const LayeredMedium &lm,
                          cmplx k_rho,
                          real z,
                          real z_,
                          EmMode mode,
                          bool direct_term = true,
                          RiemannSheet sheet = RiemannSheet::I);

// Return a SommerfeldIntegral object which can be evaluated in spatial domain.
SommerfeldIntegral get_sommerfeld_integral(const LayeredMedium &lm,
                                           real nu,
                                           EmMode mode,
                                           bool direct_term,
                                           SiParams si_params = SiParams());

} // namespace layered_media

} // namespace mthesis::gf::scalar

#endif // MTHESIS_GF_SCALAR_HPP
