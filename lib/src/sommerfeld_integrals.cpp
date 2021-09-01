#include <mthesis/detail/sommerfeld_integrals.hpp>
#include <mthesis/detail/partition_extrapolation.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cmath>
#include <cassert>

namespace mthesis
{
    cmplx si_head(const SpectralGF &gf, real nu, real rho, real a)
    {
        real w = 1.0 / rho;  // Detour width.
        cmplx z(a / 2.0, w); // Waypoint of roof-like detour.
        assert(std::isfinite(std::abs(gf.f(z))));

        auto f = [&](cmplx k_rho)
        {
            return integrand_sip(gf, nu, rho, k_rho);
        };

        using boost::math::quadrature::gauss_kronrod;
        unsigned max_depth = 15;

        auto f_1 = [&](real t)
        {
            return f(z * t);
        };
        auto val_1 = gauss_kronrod<real, 15>::integrate(f_1, 0, 1, max_depth);

        auto f_2 = [&](real t)
        {
            return f(z * (1 - t));
        };
        auto val_2 = gauss_kronrod<real, 15>::integrate(f_2, 0, 1, max_depth);

        return val_1 + val_2;
    }

    cmplx si_tail(const SpectralGF &gf,
                  real nu,
                  real rho,
                  real a,
                  PEParams params)
    {
        if (std::isfinite(gf.alpha) && std::isfinite(gf.zeta))
        {
            return pe_mosig_michalski(gf, nu, rho, a, params);
        }
        else
        {
            return pe_levin_sidi(gf, nu, rho, a, params);
        }
    }

    cmplx si_sip(const SpectralGF &gf,
                 real nu,
                 real rho,
                 real a,
                 PEParams pe_params)
    {
        if (rho == 0.0 && nu != 0.0)
        {
            return 0.0;
        }

        auto head = si_head(gf, nu, rho, a);
        auto tail = si_tail(gf, nu, rho, a, pe_params);

        return head + tail;
    }

} // namespace mthesis
