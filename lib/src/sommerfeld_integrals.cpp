#include <mthesis/sommerfeld_integrals.hpp>
#include <mthesis/partition_extrapolation.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <complex_bessel.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cmath>
#include <cassert>

namespace mthesis {

SommerfeldIntegral::SommerfeldIntegral(SpectralGF f, real nu,
                                       const LayeredMedium &lm)
    : f(f),
      nu(nu),
      lm(lm)
{}

cmplx SommerfeldIntegral::eval_integrand_sip(real rho, real z, real z_,
                                             real k_rho) const
{
    return f(z, z_, k_rho) * k_rho * boost::math::cyl_bessel_j(nu, rho * k_rho);
}

cmplx SommerfeldIntegral::eval_integrand_sip(real rho, real z, real z_,
                                             cmplx k_rho) const
{
    return f(z, z_, k_rho) * k_rho * sp_bessel::besselJ(nu, rho * k_rho);
}

cmplx SommerfeldIntegral::eval_si_along_sip(real rho, real z, real z_,
                                            pe::Params pe_params) const
{
    if (rho == 0.0 && nu != 0.0) {
        return 0.0;
    }

    auto a = calc_pe_start();

    auto head = eval_head_ellipsis(rho, z, z_, a);
    auto tail = eval_tail_on_sip(rho, z, z_, a, pe_params);

    return (head + tail) / (2.0 * M_PI);
}

cmplx SommerfeldIntegral::eval_si_along_sip(const VectorR3 &r,
                                            const VectorR3 &r_,
                                            pe::Params pe_params) const
{
    auto coords = LayeredMediumCoords(r, r_);

    return eval_si_along_sip(coords.rho, coords.z, coords.z_, pe_params);
}

real SommerfeldIntegral::calc_pe_start() const
{
    std::vector<real> k_re(lm.media.size());
    for (size_t n = 0; n < lm.media.size(); n++) {
        k_re[n] = lm.media[n].k.real();
    }

    real k_re_max = *std::max_element(k_re.begin(), k_re.end());

    real a = 1.2 * k_re_max;

    return a;
}

real SommerfeldIntegral::calc_indention(real rho, real a) const
{
    real w;
    if (rho > 0.0) {
        w = 1.0 / rho;
    } else {
        w = a / 2.0;
    }

    return w;
}

cmplx SommerfeldIntegral::eval_head_ellipsis(real rho,
                                             real z,
                                             real z_,
                                             real a) const
{
    using std::complex_literals::operator""i;

    real r = a / 2.0;
    real w = calc_indention(rho, a);

    auto f_along_path = [&](real t)
    {
        // For 0 <= t <= 1.
        real arg = M_PI * (1.0 - t);
        auto sin_val = std::sin(arg);
        auto cos_val = std::cos(arg);
        cmplx gamma = r * (1.0 + cos_val) + 1.0i * w * sin_val;
        cmplx gamma_ = r * M_PI * sin_val - 1.0i * w * M_PI * cos_val;
        return eval_integrand_sip(rho, z, z_, gamma) * gamma_;
    };

    constexpr boost::math::quadrature::gauss_kronrod<real, 31> quad;
    auto val = quad.integrate(f_along_path, 0.0, 1.0);

    return val;
}

cmplx SommerfeldIntegral::eval_tail_on_sip(real rho,
                                           real z,
                                           real z_,
                                           real a,
                                           const pe::Params &pe_params) const
{
    auto integrand = [=](real k_rho)
    {
        return eval_integrand_sip(rho, z, z_, k_rho);
    };

    cmplx val;
    if (rho > 0.0) {
        return pe::levin_sidi(integrand, nu, rho, a, pe_params);
    } else if (rho == 0 && std::abs(z - z_) > 0.0) {
        constexpr real b = std::numeric_limits<real>::infinity();
        constexpr boost::math::quadrature::gauss_kronrod<real, 15> quad;
        val = quad.integrate(integrand, a, b);
    } else {
        std::cerr << "SI with rho = |z - z_| = 0. Returning NaN.";
        val = std::numeric_limits<real>::quiet_NaN();
    }

    return val;
}

} // namespace mthesis
