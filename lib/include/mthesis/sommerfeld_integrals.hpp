#ifndef MTHESSI_SI_SOMMERFELD_INTEGRALS_HPP
#define MTHESSI_SI_SOMMERFELD_INTEGRALS_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/spectral_gf.hpp>
#include <mthesis/pe_params.hpp>

#include <functional>

namespace mthesis::si
{
    real calc_indention(const SpectralGF &gf, real rho, real a);

    cmplx eval_head_elliplis(const SpectralGF &gf, real nu, real rho, real a);

    cmplx eval_tail(const SpectralGF &gf,
                    real nu,
                    real rho,
                    real a,
                    pe::Params params);

    real get_a(const LayeredMedium &lm);

    cmplx eval_along_sip(const SpectralGF &gf,
                         real nu,
                         real rho,
                         pe::Params pe_params = pe::Params());

} // namespace mthesis::si

#endif
