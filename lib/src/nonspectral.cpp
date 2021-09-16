#include <mthesis/nonspectral.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <complex_bessel.h>
#include "./../../submodules/faddeeva/Faddeeva.hh"

#include <cassert>

namespace mthesis::nonspectral {

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

cmplx get_residue(std::function<cmplx (cmplx)> f, cmplx k_p)
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

cmplx fun_C8(cmplx p)
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

    cmplx I_q = sqrt(M_PI / rho) / 2.0 * fun_C8(p) ;

    return I_q;
}

cmplx calc_B_p(const LayeredMedium &lm,
               real n,
               std::function<cmplx (cmplx)> f,
               cmplx k_p,
               real rho)
{
    using std::complex_literals::operator""i;
    auto R_p = get_residue(f, k_p);
    auto s_p = calc_s_p(lm, k_p);
    cmplx B_p = 1.0i * R_p / (2.0 * s_p) * sqrt(k_p) *
            normalized_hankel_fun(n, k_p * rho);
    return B_p;
}

} // namespace mthesis::nonspectral
