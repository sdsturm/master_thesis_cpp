#ifndef MTHESIS_SGF_HPP
#define MTHESIS_SGF_HPP

#include <mthesis/detail/definitions.hpp>
#include <mthesis/detail/solution_domain.hpp>
#include <mthesis/detail/tlgf.hpp>
#include <mthesis/detail/spectral_gf.hpp>

// Generic scalar Green's functions.
namespace mthesis::sgf
{
    cmplx free_space(const Medium &m, const VectorR3 &r, const VectorR3 &r_);

    cmplx lm_generic_spectral(const LayeredMedium &lm,
                              real z,
                              real z_,
                              cmplx k_rho,
                              EmMode mode,
                              bool direct_term);

    const si::SpectralGF lm_get_generic_spectral_gf(const LayeredMedium &lm,
                                                    const VectorR3 &r,
                                                    const VectorR3 &r_,
                                                    EmMode mode,
                                                    bool direct_term);

    cmplx lm_generic_spatial(const LayeredMedium &lm,
                             const VectorR3 &r,
                             const VectorR3 &r_,
                             EmMode mode,
                             bool direct_term);

} // namespace mthesis::sgf

#endif
