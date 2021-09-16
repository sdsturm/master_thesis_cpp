#include <mthesis/nonspectral.hpp>

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

} // namespace mthesis::nonspectral
