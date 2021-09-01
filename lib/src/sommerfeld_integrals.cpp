#include <mthesis/detail/sommerfeld_integrals.hpp>
#include <mthesis/detail/partition_extrapolation.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cmath>
#include <cassert>

namespace mthesis::si
{
    cmplx eval_head_rooftop(const SpectralGF &gf, real nu, real rho, real a)
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

    cmplx eval_tail(const SpectralGF &gf,
                    real nu,
                    real rho,
                    real a,
                    pe::Params params)
    {
        auto f = [&](real k_rho)
        {
            return integrand_sip(gf, nu, rho, k_rho);
        };

        cmplx val;
        if (rho > 0.0)
        {
            if (std::isfinite(gf.alpha) && std::isfinite(gf.zeta))
            {
                return pe::mosig_michalski(f, gf.alpha, gf.zeta, nu, rho, a, params);
            }
            else
            {
                return pe::levin_sidi(f, nu, rho, a, params);
            }
        }
        else if (rho == 0 && std::abs(gf.z - gf.z_) > 0.0)
        {
            val = boost::math::quadrature::gauss_kronrod<real, 15>::integrate(
                f, a, std::numeric_limits<real>::infinity(), 5);
        }
        else
        {
            std::cerr << "SI with rho = |z - z_| = 0. Returning NaN.";
            val = std::numeric_limits<real>::quiet_NaN();
        }
        return val;
    }

    real get_a(const LayeredMedium &lm)
    {
        std::vector<real> k_re(lm.media.size());
        for (size_t n = 0; n < lm.media.size(); n++)
            k_re[n] = lm.media[n].k.real();

        real k_re_max = *std::max_element(k_re.begin(), k_re.end());

        real a = 1.2 * k_re_max;

        return a;
    }

    cmplx eval_along_sip(const SpectralGF &gf,
                         real nu,
                         real rho,
                         pe::Params pe_params)
    {
        if (rho == 0.0 && nu != 0.0)
        {
            return 0.0;
        }

        auto a = get_a(gf.lm);

        auto head = eval_head_rooftop(gf, nu, rho, a);
        auto tail = eval_tail(gf, nu, rho, a, pe_params);

        return (head + tail) / (2.0 * M_PI);
    }

} // namespace mthesis::si
