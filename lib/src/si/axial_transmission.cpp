#include <mthesis/si/axial_transmission.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace mthesis::si::axial_transmission {

cmplx eval_si_along_sip(const SommerfeldIntegral &si,
                        real rho,
                        real z,
                        real z_,
                        pe::Params pe_params)
{
    if (rho == 0.0 && si.nu != 0.0) {
        return 0.0;
    }

    auto a = utils::calc_pe_start(si);

    auto head = utils::eval_head_ellipsis(si, rho, z, z_, a);
    auto tail = utils::eval_tail_on_sip(si, rho, z, z_, a, pe_params);

    return (head + tail) / (2.0 * M_PI);
}

cmplx eval_si_along_sip(const SommerfeldIntegral &si,
                        const VectorR3 &r,
                        const VectorR3 &r_,
                        pe::Params pe_params)
{
    LMCoords c(r, r_);
    return eval_si_along_sip(si, c.rho, c.z, c.z_, pe_params);
}

namespace utils {

real calc_pe_start(const SommerfeldIntegral &si)
{
    real k_re_max = 0;
    for (const auto &medium : si.lm.media) {
        if (medium.k.real() > k_re_max) {
            k_re_max = medium.k.real();
        }
    }

    // Note: something like 3 * k_re_max slows down the algorithm significantly.
    real a = 1.2 * k_re_max;  // Empirical.

    return a;
}

real calc_indention(real rho, real a)
{
    real w;
    if (rho > 0.0) {
        w = 1.0 / rho;
    } else {
        w = a / 2.0;
    }

    return w;
}

cmplx eval_head_ellipsis(const SommerfeldIntegral &si,
                         real rho,
                         real z,
                         real z_,
                         real a)
{
    using std::complex_literals::operator""i;

    real r = a / 2.0;
    real w = calc_indention(rho, a);

    auto f = [=](cmplx k_rho)
    {
        return si.eval_integrand_sip(k_rho, rho, z, z_, RiemannSheet::I);
    };

    auto f_along_path = [&](real t)
    {
        // For 0 <= t <= 1.
        real arg = M_PI * (1.0 - t);
        auto sin_val = std::sin(arg);
        auto cos_val = std::cos(arg);
        cmplx gamma = r * (1.0 + cos_val) + 1.0i * w * sin_val;
        cmplx gamma_ = r * M_PI * sin_val - 1.0i * w * M_PI * cos_val;
        return f(gamma) * gamma_;
    };

    constexpr boost::math::quadrature::gauss_kronrod<real, 31> quad;
    auto val = quad.integrate(f_along_path, 0.0, 1.0);

    return val;
}

cmplx eval_tail_on_sip(const SommerfeldIntegral &si,
                       real rho,
                       real z,
                       real z_,
                       real a,
                       const pe::Params &pe_params)
{
    auto integrand = [=](real k_rho)
    {
        return si.eval_integrand_sip(k_rho, rho, z, z_, RiemannSheet::I);
    };

    cmplx val;
    if (rho > 0.0) {
        return pe::levin_sidi(integrand, si.nu, rho, a, pe_params);
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

}

// namespace utils

} // namespace mthesis::si::axial_transmission
