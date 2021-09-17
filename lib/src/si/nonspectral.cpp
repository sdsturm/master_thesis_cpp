#include <mthesis/si/nonspectral.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <complex_bessel.h>
#include "./../../submodules/faddeeva/Faddeeva.hh"

#include <cassert>

namespace mthesis::si::nonspectral {

// Hide implementation details in nested namespace.
namespace utils {

cmplx normalized_hankel(real n, cmplx z)
{
    using std::complex_literals::operator""i;
    return pow(-1.0i, n) * sqrt(M_PI * z / 2.0i) *
            exp(1.0i * z) * sp_bessel::sph_hankelH2(n, z);
}

// Hide implementation details in nested namespace.
namespace tm {

cmplx get_k_p(const LayeredMedium &lm)
{
    assert(lm.is_sommerfeld_half_space());

    cmplx k_p = lm.media.back().k * sqrt(lm.media.front().eps_r /
                                             (1.0 + lm.media.front().eps_r));
    return k_p;
}

cmplx get_s_p(const LayeredMedium &lm, cmplx k_p)
{
    using std::complex_literals::operator""i;

    cmplx s_p = sqrt(1.0i * (k_p - lm.media.back().k));

    return s_p;
}

cmplx calc_k_p_residue(const SommerfeldIntegral &si, real z, real z_, cmplx k_p)
{
    using std::complex_literals::operator""i;
    constexpr boost::math::quadrature::gauss_kronrod<real, 31> quad;

    auto f = [=](real t)
    {
    // Curve and derivative for 0 <= t <= 1.
    real r = 0.5 * abs(k_p.imag()); 	// Just a try.
    cmplx gamma =  k_p + r * exp(2.0i * M_PI * t);
    cmplx gamma_ =  2.0i * M_PI * r * exp(2.0i * M_PI * t);
    return si.eval_spectral_gf(z, z_, gamma) * gamma_;
    };

    cmplx R_p = quad.integrate(f, 0.0, 1.0);
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

cmplx integrand_I_p(cmplx s, cmplx s_p, cmplx B_p, real rho)
{
    // TODO
    return 0;
}

} // namespace tm

// Hide implementation details in nested namespace.
namespace te {

}

// namespace te

} // namespace utils

} // namespace mthesis::si::nonspectral
