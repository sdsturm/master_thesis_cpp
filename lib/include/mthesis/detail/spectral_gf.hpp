#ifndef MTHESIS_SPECTRAL_GF_HPP
#define MTHESIS_SPECTRAL_GF_HPP

#include <mthesis/detail/definitions.hpp>
#include <mthesis/detail/solution_domain.hpp>

#include <functional>

namespace mthesis
{
    struct SpectralGF
    {
        const std::function<cmplx(cmplx)> f;
        const real alpha;
        const real zeta;
        const std::vector<cmplx> bp;  // Branch points.
        const std::vector<cmplx> swp; // Surface wave poles.

        SpectralGF(std::function<cmplx(cmplx)> f,
                   const LayeredMedium &lm,
                   real alpha = std::numeric_limits<real>::quiet_NaN(), // Flag.
                   real zeta = std::numeric_limits<real>::quiet_NaN(),  // Flag.
                   bool identify_poles = false);
    };

    std::vector<cmplx> get_branch_points(const LayeredMedium &lm);

    std::vector<cmplx> identify_swp(std::function<cmplx(cmplx)> f);

    cmplx integrand_sip(const SpectralGF &gf, real nu, real rho, real k_rho);

    cmplx integrand_sip(const SpectralGF &gf, double nu, real rho, cmplx k_rho);

    cmplx integrand_eip(const SpectralGF &gf, double nu, real rho, cmplx k_rho);

} // namespace mthesis

#endif
