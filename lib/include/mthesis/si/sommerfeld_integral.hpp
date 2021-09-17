#ifndef MTHESSI_SOMMERFELD_INTEGRAL_HPP
#define MTHESSI_SOMMERFELD_INTEGRAL_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/si/partition_extrapolation.hpp>

#include <functional>

namespace mthesis {

struct SommerfeldIntegral
{
    SommerfeldIntegral(SpectralGF f, real nu, const LayeredMedium &lm);

    cmplx eval_spectral_gf(real z, real z_, cmplx k_rho) const;
    cmplx eval_integrand_sip(real rho, real z, real z_, real k_rho) const;
    cmplx eval_integrand_sip(real rho, real z, real z_, cmplx k_rho) const;

    const SpectralGF f;
    const real nu;
    const LayeredMedium &lm;
};

} // namespace mthesis

#endif // MTHESSI_SOMMERFELD_INTEGRAL_HPP
