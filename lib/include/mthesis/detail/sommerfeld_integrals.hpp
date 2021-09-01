#ifndef MTHESSI_SOMMERFELD_INTEGRALS_HPP
#define MTHESSI_SOMMERFELD_INTEGRALS_HPP

#include <mthesis/detail/definitions.hpp>
#include <mthesis/detail/solution_domain.hpp>
#include <mthesis/detail/spectral_gf.hpp>
#include <mthesis/detail/pe_params.hpp>

#include <functional>

namespace mthesis
{
    cmplx si_head(const SpectralGF &gf, real nu, real rho, real a);

    cmplx si_tail(const SpectralGF &gf,
                  real nu,
                  real rho,
                  real a,
                  PEParams params);

    cmplx si_sip(const SpectralGF &gf,
                 real nu,
                 real rho,
                 real a,
                 PEParams pe_params);

} // namespace mthesis

#endif
