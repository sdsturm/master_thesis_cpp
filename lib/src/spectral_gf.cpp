#include <mthesis/detail/spectral_gf.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <complex_bessel.h>

namespace mthesis::si
{
    std::vector<cmplx> get_branch_points(const LayeredMedium &lm)
    {
        std::vector<cmplx> bp_locations;

        if (std::isinf(lm.d.front()))
            bp_locations.push_back(lm.media.front().k);

        if (std::isinf(lm.d.back()))
            bp_locations.push_back(lm.media.back().k);

        return bp_locations;
    }

    std::vector<cmplx> identify_swp(/* TODO */)
    {
        std::vector<cmplx> sip_locations;

        // TODO

        return sip_locations;
    }

    SpectralGF::SpectralGF(std::function<cmplx(cmplx)> f,
                           real z,
                           real z_,
                           const LayeredMedium &lm,
                           real alpha,
                           real zeta,
                           bool identify_poles)
        : f(f),
          z(z),
          z_(z_),
          lm(lm),
          alpha(alpha),
          zeta(zeta),
          bp(identify_poles ? get_branch_points(lm) : std::vector<cmplx>(0)),
          swp(identify_poles ? identify_swp() : std::vector<cmplx>(0))
    {
    }

    SpectralGF::SpectralGF(std::function<cmplx(cmplx)> f,
                           real z,
                           real z_,
                           const LayeredMedium &lm,
                           bool identify_poles)
        : SpectralGF(f,
                     z,
                     z_,
                     lm,
                     std::numeric_limits<real>::quiet_NaN(),
                     std::numeric_limits<real>::quiet_NaN(),
                     identify_poles)
    {
    }

    cmplx integrand_sip(const SpectralGF &gf, real nu, real rho, real k_rho)
    {
        return boost::math::cyl_bessel_j(nu, rho * k_rho) * k_rho * gf.f(k_rho);
    }

    cmplx integrand_sip(const SpectralGF &gf, double nu, real rho, cmplx k_rho)
    {
        return sp_bessel::besselJ(nu, rho * k_rho) * k_rho * gf.f(k_rho);
    }

    cmplx integrand_eip(const SpectralGF &gf, double nu, real rho, cmplx k_rho)
    {
        return sp_bessel::hankelH2(nu, rho * k_rho) * k_rho * gf.f(k_rho);
    }

} // namespace mthesis::si
