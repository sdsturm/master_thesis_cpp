#ifndef MTHESIS_SI_SPECTRAL_GF_HPP
#define MTHESIS_SI_SPECTRAL_GF_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

#include <functional>

namespace mthesis::si
{
    std::vector<cmplx> get_branch_points(const LayeredMedium &lm);

    std::vector<cmplx> identify_swp(/* TODO */);

    struct SpectralGF
    {
        const std::function<cmplx(cmplx)> f;
        const real z;
        const real z_;
        const LayeredMedium &lm;
        const real alpha;
        const real zeta;
        const std::vector<cmplx> bp;  // Branch points.
        const std::vector<cmplx> swp; // Surface wave poles.

        SpectralGF(std::function<cmplx(cmplx)> f,
                   real z,
                   real z_,
                   const LayeredMedium &lm,
                   real alpha, // NaN if not given.
                   real zeta,  // NaN if not given.
                   bool identify_poles = false);

        SpectralGF(std::function<cmplx(cmplx)> f,
                   real z,
                   real z_,
                   const LayeredMedium &lm,
                   bool identify_poles = false);
    };

    cmplx integrand_sip(const SpectralGF &gf, real nu, real rho, real k_rho);

    cmplx integrand_sip(const SpectralGF &gf, double nu, real rho, cmplx k_rho);

} // namespace mthesis::si

#endif
