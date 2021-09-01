#ifndef MTHESSI_SI_SOMMERFELD_INTEGRALS_HPP
#define MTHESSI_SI_SOMMERFELD_INTEGRALS_HPP

#include <mthesis/detail/definitions.hpp>
#include <mthesis/detail/solution_domain.hpp>
#include <mthesis/detail/spectral_gf.hpp>
#include <mthesis/detail/pe_params.hpp>

#include <functional>

namespace mthesis::si
{
    cmplx eval_head_rooftop(const SpectralGF &gf, real nu, real rho, real a);

    cmplx eval_tail(const SpectralGF &gf,
                    real nu,
                    real rho,
                    real a,
                    pe::Params params);

    cmplx eval_along_sip(const SpectralGF &gf,
                         real nu,
                         real rho,
                         real a,
                         pe::Params params);

} // namespace mthesis::si

#endif
