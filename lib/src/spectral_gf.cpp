#include <mthesis/detail/spectral_gf.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <complex_bessel.h>

namespace mthesis
{
    std::vector<cmplx> get_branch_points(const LayeredMedium &lm)
    {
        std::vector<cmplx> bp_locations;

        return bp_locations;
    }

    std::vector<cmplx> identify_swp(std::function<cmplx(cmplx)> f)
    {
        std::vector<cmplx> sip_locations;

        return sip_locations;
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

    SpectralGF::SpectralGF(std::function<cmplx(cmplx)> f,
                           const LayeredMedium &lm,
                           real alpha,
                           real zeta,
                           bool identify_poles)
        : f(f),
          alpha(alpha),
          zeta(zeta),
          bp(identify_poles ? get_branch_points(lm) : std::vector<cmplx>(0)),
          swp(identify_poles ? identify_swp(f) : std::vector<cmplx>(0))
    {
    }
    
} // namespace mthesis
