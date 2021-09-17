#include <mthesis/nonspectral.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <complex_bessel.h>
#include "./../../submodules/faddeeva/Faddeeva.hh"

#include <cassert>

namespace mthesis::nonspectral {

#if 0
cmplx get_sommerfeld_pole(const LayeredMedium &lm, EmMode mode)
{
    // Only Sommerfeld pole in case of Sommerfeld half-space is implemented.
    // See section 2.3 in Michalski2016b.
    assert(mode == EmMode::TM);

    assert(lm.media.size() == 2);
    assert(lm.top_bc == BC::open);
    assert(lm.bottom_bc == BC::open);

    cmplx k_p = lm.media.back().k * sqrt(lm.media.front().eps_r /
                                         (1.0 + lm.media.front().eps_r));
    return k_p;
}

cmplx calc_s_p(const LayeredMedium &lm, cmplx k_p)
{
    using std::complex_literals::operator""i;
    return sqrt(1.0i * (k_p - lm.media.back().k));
}

cmplx integrate_k_p_residue(std::function<cmplx (cmplx)> f, cmplx k_p)
{
    constexpr boost::math::quadrature::gauss_kronrod<real, 31> quad;

    // Curve and derivative for 0 <= t <= 1.
    real r = 0.5 * abs(k_p.imag()); 	// Just a try.
    auto gamma = [=](real t) { return k_p + r * exp(2.0i * M_PI * t); };
    auto gamma_ = [=](real t) { return 2.0i * M_PI * r * exp(2.0i * M_PI * t); };

    auto integrand = [=](real t) { return f(gamma(t)) * gamma_(t); };
    cmplx R_p = quad.integrate(integrand, 0.0, 1.0);
    return R_p;
}

cmplx normalized_hankel_fun(real n, cmplx z)
{
    using std::complex_literals::operator""i;
    return pow(1.0i, n) * sqrt(M_PI * z / 2.0i) *
            exp(1.0i * z) * sp_bessel::sph_hankelH2(n, z);
}

cmplx cmplx_calc_num_dist_eqC5(const LayeredMedium &lm, cmplx k_p, real rho)
{
    auto s_p = calc_s_p(lm, k_p);
    return pow(s_p, 2) * rho;
}

cmplx calc_F_eqC8(cmplx p)
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

cmplx calc_I_q(const LayeredMedium &lm, cmplx k_p, real rho)
{
    using std::complex_literals::operator""i;

    cmplx p = 1.0i * (k_p - lm.media.back().k) * rho;

    cmplx I_q = sqrt(M_PI / rho) / 2.0 * calc_F_eqC8(p) ;

    return I_q;
}

cmplx calc_B_p(const LayeredMedium &lm,
               real n,
               std::function<cmplx (cmplx)> f,
               cmplx k_p,
               real rho)
{
    using std::complex_literals::operator""i;
    auto R_p = integrate_k_p_residue(f, k_p);
    auto s_p = calc_s_p(lm, k_p);
    cmplx B_p = 1.0i * R_p / (2.0 * s_p) * sqrt(k_p) *
            normalized_hankel_fun(n, k_p * rho);
    return B_p;
}

cmplx calc_f_folded(cmplx k_rho,
                    const std::function<cmplx (cmplx, RiemannSheet)> &f,
                    const LayeredMedium &lm)
{
    const cmplx &k_N = lm.media.front().k;

    // Vertical branch cut for top (air) half-space.
    assert(k_rho.real() == lm.media.back().k.real());

    real bc_k_rho_imag = k_N.real() * k_N.imag() / lm.media.back().k.real();

    RiemannSheet sheet_f_left, sheet_f_right;
    if (k_rho.imag() >= bc_k_rho_imag) {
        sheet_f_left = RiemannSheet::II;
    } else {
        sheet_f_left = RiemannSheet::IV;
    }

    if (k_rho.imag() >= bc_k_rho_imag) {
        sheet_f_right = RiemannSheet::I;
    } else {
        sheet_f_right = RiemannSheet::III;
    }

    cmplx f_left = f(k_rho, sheet_f_left);
    cmplx f_right = f(k_rho, sheet_f_right);

    return f_left - f_right;
}

cmplx calc_F_eq81(cmplx s,
                  int n,
                  const LayeredMedium &lm,
                  const std::function<cmplx (cmplx, RiemannSheet)> &f,
                  real rho)
{
    using std::complex_literals::operator""i;
    cmplx k_rho = lm.media.back().k - 1.0i * pow(s, 2);
    cmplx f_folded = calc_f_folded(k_rho, f, lm);
    cmplx H_normalized = normalized_hankel_fun(n, k_rho * rho);

    return f_folded * sqrt(k_rho) * H_normalized;
}

cmplx calc_I_p(int n,
               const LayeredMedium &lm,
               const std::function<cmplx (cmplx, RiemannSheet)> &f,
               real rho,
               cmplx s_p,
               cmplx B_p)
{
    auto integrand = [=](real s)
    {
        cmplx t1 = calc_F_eq81(s, n, lm, f, rho);
        cmplx t2 = 2.0 * B_p * s / (pow(s, 2) - pow(s_p, 2));
        cmplx t3 = exp(-pow(s, 2) * rho) * s;
        return (t1 + t2) * t3;
    };

    constexpr boost::math::quadrature::gauss_kronrod<real, 61> quad;

    return quad.integrate(integrand, 0, std::numeric_limits<real>::infinity());
}

cmplx calc_I_1_n(const LayeredMedium &lm,
                 int n,
                 const std::function<cmplx (cmplx, RiemannSheet)> &f,
                 real rho)
{
    using std::complex_literals::operator""i;
    cmplx t1 = 2.0 * pow(1.0i, n + 1) * sqrt(2.0i) *
            exp(-1.0i * lm.media.back().k * rho) / rho;

    cmplx s_p = calc_s_p(lm, k_p);

    cmplx t2 = sqrt(rho / M_PI) * calc_I_p(n, lm, f, rho, s_p, B_p);

    return val;
}

cmplx eval_nonspectral(const LayeredMedium &lm,
                       const std::function<cmplx (cmplx, RiemannSheet)> &f,
                       real rho)
{
    assert(rho > 0);
    assert(lm.is_sommerfeld_half_space());


    return 0;
}
#endif

} // namespace mthesis::nonspectral
