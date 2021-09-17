#ifndef MTHESSI_SI_SOMMERFELD_INTEGRALS_HPP
#define MTHESSI_SI_SOMMERFELD_INTEGRALS_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/partition_extrapolation.hpp>

#include <functional>

namespace mthesis {

class SommerfeldIntegral
{
public:
    using spectral_gf = std::function<cmplx(real z, real z_, cmplx k_rho)>;

    SommerfeldIntegral(spectral_gf f, real nu, const LayeredMedium &lm);

    cmplx eval_si_along_sip(real rho, real z, real z_,
                            pe::Params pe_params = pe::Params()) const;

    cmplx eval_si_along_sip(const VectorR3 &r, const VectorR3 &r_,
                            pe::Params pe_params = pe::Params()) const;

private:
    cmplx eval_integrand_sip(real rho, real z, real z_, real k_rho) const;
    cmplx eval_integrand_sip(real rho, real z, real z_, cmplx k_rho) const;

    real calc_pe_start() const;

    real calc_indention(real rho, real a) const;

    cmplx eval_head_ellipsis(real rho, real z, real z_, real a) const;

    cmplx eval_tail_on_sip(real rho, real z, real z_,
                           real a, const pe::Params &pe_params) const;

public:
    const spectral_gf f;
    const real nu;
    const LayeredMedium &lm;
};

} // namespace mthesis

#endif
