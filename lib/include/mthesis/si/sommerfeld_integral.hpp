#ifndef MTHESSI_SOMMERFELD_INTEGRAL_HPP
#define MTHESSI_SOMMERFELD_INTEGRAL_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

#include <functional>

namespace mthesis {

struct SommerfeldIntegral
{
    using SpectralGF =
    std::function<cmplx (cmplx k_rho, real z, real z_, RiemannSheet sheet)>;

    SommerfeldIntegral(SpectralGF f, real nu, const LayeredMedium &lm);

    cmplx eval_spectral_gf(cmplx k_rho,
                           real z,
                          real z_,
                           RiemannSheet sheet) const;

    cmplx eval_integrand_sip(real k_rho,
                             real rho,
                             real z,
                             real z_,
                             RiemannSheet sheet) const;

    cmplx eval_integrand_sip(cmplx k_rho,
                             real rho,
                             real z,
                             real z_,
                             RiemannSheet sheet) const;

    const SpectralGF f;
    const real nu;
    const LayeredMedium &lm;
};

} // namespace mthesis

#endif // MTHESSI_SOMMERFELD_INTEGRAL_HPP
