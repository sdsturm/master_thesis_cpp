#ifndef MTHESSI_SI_SOMMERFELD_INTEGRALS_HPP
#define MTHESSI_SI_SOMMERFELD_INTEGRALS_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/partition_extrapolation.hpp>

#include <functional>

namespace mthesis {

struct SiParams
{
    real alpha;				// Michalski2016a (70)
    real zeta;				// Michalski2016a (70)
    bool identify_singularities;

    SiParams(real alpha, real zeta, bool identify_singularities = false);
    SiParams();
};


class SommerfeldIntegral
{
public:
    using spectral_gf = std::function<cmplx(real z, real z_, cmplx k_rho)>;

    SommerfeldIntegral(spectral_gf f,
                       real nu,
                       const LayeredMedium &lm,
                       SiParams params = SiParams());

    cmplx eval_spectral_gf(real z, real z_, cmplx k_rho) const;

    // Note: prefactor 1 / (2 * pi) for SIP.
    cmplx eval_integrand_sip(real rho, real z, real z_, real k_rho) const;
    cmplx eval_integrand_sip(real rho, real z, real z_, cmplx k_rho) const;

    cmplx eval_si_along_sip(real rho, real z, real z_,
                            pe::Params pe_params = pe::Params()) const;

    cmplx eval_si_along_sip(const VectorR3 &r, const VectorR3 &r_,
                            pe::Params pe_params = pe::Params()) const;

    // Note: prefactor 1 / (4 * pi) for EIP.
    cmplx eval_integrand_eip(real rho, real z, real z_, cmplx k_rho) const;

    const LayeredMedium &get_lm() const;

private:
    real calc_pe_start() const;

    real calc_indention(real rho, real a) const;

    cmplx eval_head_ellipsis(real rho, real z, real z_, real a) const;

    cmplx eval_tail_on_sip(real rho, real z, real z_,
                           real a, const pe::Params &pe_params) const;

private:
    spectral_gf f;
    real nu;
    const LayeredMedium &lm;
    const SiParams params;
    std::vector<cmplx> bp; 	// Branch points
    std::vector<cmplx> swp;	// Surface wave poles
};

std::vector<cmplx> get_branch_points(const LayeredMedium &lm);

std::vector<cmplx> identify_swp(/* TODO */);

} // namespace mthesis

#endif
