#include <mthesis/si/nonspectral.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <complex_bessel.h>
#include "./../../submodules/faddeeva/Faddeeva.hh"

#include <cassert>

namespace mthesis::si::nonspectral {

cmplx eval_nonspectral(const SommerfeldIntegral &si,
                       const LMCoords &c,
                       EmMode mode)
{
    cmplx val;

    switch (mode) {
    case EmMode::TM:
        val =  utils::tm::calc_I1_eq85(si, c);
        break;
    case EmMode::TE:
        // TODO
        break;
    }

    return val;
}

// =============================================================================

namespace utils {

cmplx normalized_hankel(real n, cmplx z)
{
    using std::complex_literals::operator""i;
    return pow(-1.0i, n) * sqrt(M_PI * z / 2.0i) *
            exp(1.0i * z) * sp_bessel::sph_hankelH2(n, z);
}

cmplx calc_folded_sgf(const SommerfeldIntegral &si,
                      const LMCoords &c,
                      cmplx k_rho)
{
    const cmplx &k_1 = si.lm.media.back().k;
    const cmplx &k_N = si.lm.media.front().k;
    real k_rho_im_intersect = k_N.real() * k_N.imag() / k_1.real();

    RiemannSheet sheet_left, sheet_right;
    if (k_rho.imag() > k_rho_im_intersect) {
        sheet_left = RiemannSheet::II;
        sheet_right = RiemannSheet::IV;
    } else {
        sheet_left = RiemannSheet::I;
        sheet_right = RiemannSheet::III;
    }

    // TODO: incorporate Riemann sheet.
    cmplx f_left = si.eval_spectral_gf(k_rho, c.z, c.z_, sheet_left);
    cmplx f_right = si.eval_spectral_gf(k_rho, c.z, c.z_, sheet_right);

    return f_left - f_right;
}

cmplx calc_F_eq81(cmplx s, const SommerfeldIntegral &si, const LMCoords &c)
{
    using std::complex_literals::operator""i;

    cmplx k_rho = si.lm.media.back().k - 1.0i * pow(s, 2);

    return calc_folded_sgf(si, c, k_rho) * sqrt(k_rho) *
            normalized_hankel(si.nu, k_rho * c.rho);
}

// =============================================================================

namespace tm {

cmplx calc_k_p(const LayeredMedium &lm)
{
    assert(lm.is_sommerfeld_half_space());
    assert(lm.media.back().k.real() <= lm.media.front().k.real());

    cmplx k_p = lm.media.back().k * sqrt(lm.media.front().eps_r /
                                             (1.0 + lm.media.front().eps_r));
    return k_p;
}

cmplx calc_s_p(const LayeredMedium &lm, cmplx k_p)
{
    using std::complex_literals::operator""i;

    cmplx s_p = sqrt(1.0i * (k_p - lm.media.back().k));

    return s_p;
}

#if 0
cmplx calc_R_p_quad(const SommerfeldIntegral &si, const LMCoords &c, cmplx k_p)
{
    using std::complex_literals::operator""i;
    constexpr boost::math::quadrature::gauss_kronrod<real, 31> quad;

    auto f = [=](real t)
    {
    // Curve and derivative for 0 <= t <= 1.
    real r = 0.5 * abs(k_p.imag()); 	// Just a try.
    cmplx gamma =  k_p + r * exp(2.0i * M_PI * t);
    cmplx gamma_ =  2.0i * M_PI * r * exp(2.0i * M_PI * t);
    return si.eval_spectral_gf(c.z, c.z_, gamma) * gamma_;
    };

    cmplx R_p = quad.integrate(f, 0.0, 1.0);
    return R_p;
}
#endif

cmplx calc_R_p_closed_form(const LayeredMedium &lm, const LMCoords &c)
{
    using std::complex_literals::operator""i;

    const cmplx &eps_r = lm.media.front().eps_r;
    const cmplx &k_1 = lm.media.back().k;

    cmplx R_p = -1.0i * eps_r * sqrt(eps_r) / (pow(eps_r, 2) - 1.0) *
            exp(1.0i * k_1 * abs(c.z - c.z_) / (eps_r + 1.0));

    return R_p;
}

cmplx calc_numerical_distance(const LayeredMedium &lm, cmplx k_p, real rho)
{
    using std::complex_literals::operator""i;

    cmplx p = 1.0i * (k_p - lm.media.back().k) * rho;

    return p;
}

cmplx calc_F_C8(cmplx p)
{
    using std::complex_literals::operator""i;

    cmplx imagterm = 1.0i * sqrt(M_PI * p);
    if (p.imag() > 0.0) {
        imagterm *= Faddeeva::w(sqrt(p));
    } else if (p.imag() < 0.0) {
        imagterm *= Faddeeva::w(sqrt(p)) - exp(-p);
    } else {
        imagterm *= -Faddeeva::w(-sqrt(p));
    }

    return 1.0 + imagterm;
}

cmplx calc_B_p(cmplx k_p, cmplx s_p, cmplx R_p, real n, real rho)
{
    using std::complex_literals::operator""i;

    cmplx B_p = 1.0i * R_p / (2.0 * s_p) * sqrt(k_p) *
            normalized_hankel(n, k_p * rho);

    return B_p;
}

cmplx integrand_I_p(const SommerfeldIntegral &si, const LMCoords &c, cmplx s,
                    cmplx s_p, cmplx B_p)
{
    cmplx t1 = calc_F_eq81(s, si, c);
    cmplx t2 = 2.0 * B_p * s / (pow(s, 2) - pow(s_p, 2));
    cmplx t3 = exp(-pow(s, 2) * c.rho) * s;

    return (t1 + t2) * t3;
}

cmplx calc_I_p(const SommerfeldIntegral &si, const LMCoords &c,
               cmplx B_p, cmplx s_p)
{
    auto f = [=](real s)
    {
        return integrand_I_p(si, c, s, s_p, B_p);
    };

    constexpr boost::math::quadrature::gauss_kronrod<real, 61> quad;

    auto val = quad.integrate(f, 0, std::numeric_limits<real>::infinity());

    return val;
}

cmplx calc_I1_eq85(const SommerfeldIntegral &si, const LMCoords &c)
{
    using std::complex_literals::operator""i;

    const cmplx &k_1 = si.lm.media.back().k;
    cmplx k_p = calc_k_p(si.lm);
    cmplx s_p = calc_s_p(si.lm, k_p);
    cmplx R_p= calc_R_p_closed_form(si.lm, c);
//    cmplx R_p = calc_R_p_quad(si, c, k_p);

    cmplx p = calc_numerical_distance(si.lm, k_p, c.rho);
    cmplx B_p = calc_B_p(k_p, s_p, R_p, si.nu, c.rho);

    cmplx I_p = calc_I_p(si, c, B_p, s_p);

    cmplx val = 2.0 * pow(1.0i, si.nu + 1.0) * sqrt(2.0i) *
            (sqrt(c.rho / M_PI) * I_p - B_p * calc_F_C8(p)) *
            exp(-1.0i * k_1 * c.rho) / c.rho;

    return val;
}

} // namespace tm

// =============================================================================

namespace te {

cmplx integrand_I_p(const SommerfeldIntegral &si, const LMCoords &c, cmplx s)
{
    cmplx t1 = calc_F_eq81(s, si, c);
    cmplx t3 = exp(-pow(s, 2) * c.rho) * s;

    return t1 * t3;
}

cmplx calc_I_p(const SommerfeldIntegral &si, const LMCoords &c)
{
    auto f = [=](real s)
    {
        return integrand_I_p(si, c, s);
    };

    constexpr boost::math::quadrature::gauss_kronrod<real, 61> quad;

    return quad.integrate(f, 0, std::numeric_limits<real>::infinity());
}

cmplx calc_I1_eq85(const SommerfeldIntegral &si, const LMCoords &c)
{
    using std::complex_literals::operator""i;

    const cmplx &k_1 = si.lm.media.back().k;

    cmplx I_p = calc_I_p(si, c);

    cmplx val = 2.0 * pow(1.0i, si.nu + 1.0) * sqrt(2.0i) *
            (sqrt(c.rho / M_PI) * I_p) *
            exp(-1.0i * k_1 * c.rho) / c.rho;

    return val;
}

} // namespace te

} // namespace utils

} // namespace mthesis::si::nonspectral
