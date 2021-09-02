#include <mthesis/sommerfeld_integrals.hpp>
#include <mthesis/partition_extrapolation.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cmath>
#include <cassert>

namespace mthesis::si
{
    real calc_indention(const SpectralGF &gf, real rho, real a)
    {
        real w;
        if (rho > 0.0)
            w = 1.0 / rho;
        else
            w = a / 2.0;

        auto test_val = gf.f(cmplx(a / 2.0, w));
        assert(std::isfinite(std::abs(test_val)));

        return w;
    }

    cmplx eval_head_elliplis(const SpectralGF &gf, real nu, real rho, real a)
    {
        using std::complex_literals::operator""i;

        auto f = [&](cmplx k_rho)
        {
            return integrand_sip(gf, nu, rho, k_rho);
        };

        real r = a / 2.0;
        real w = calc_indention(gf, rho, a);

        auto f_e = [&](real t)
        {
            real arg = M_PI * (1.0 - t);
            auto sin_val = std::sin(arg);
            auto cos_val = std::cos(arg);
            cmplx gamma = r * (1.0 + cos_val) + 1.0i * w * sin_val;
            cmplx gamma_ = r * M_PI * sin_val - 1.0i * w * M_PI * cos_val;
            return f(gamma) * gamma_;
        };

        constexpr boost::math::quadrature::gauss_kronrod<real, 31> quad;
        auto val = quad.integrate(f_e, 0.0, 1.0);

        return val;
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
            constexpr boost::math::quadrature::gauss_kronrod<real, 15> quad;
            val = quad.integrate(f, a, std::numeric_limits<real>::infinity());
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

        auto head = eval_head_elliplis(gf, nu, rho, a);
        auto tail = eval_tail(gf, nu, rho, a, pe_params);

        return (head + tail) / (2.0 * M_PI);
    }

} // namespace mthesis::si
