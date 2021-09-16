#include <mthesis/nonspectral.hpp>

#include <complex_bessel.h>

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

cmplx faddeeva(cmplx z)
{
    return exp(-pow(z, 2));
}

cmplx func_C8(cmplx p)
{
    using std::complex_literals::operator""i;
    cmplx val = 1.0 + 1.0i * sqrt(M_PI * p);
    if (p.imag() > 0.0) {
        val *= 0;
    } else if (p.imag() < 0.0) {
        val *= 0;
    } else {
        val *= 0;
    }
    return val;
}

#if 0
cmplx calc_I_q(const LayeredMedium &lm, cmplx k_p, real rho)
{
    using std::complex_literals::operator""i;
//    cmplx p = 1.0i * (k_p - lm.media.back().k) * rho;

    cmplx I_q = 0.5 * sqrt(M_PI / rho);
    return I_q;
}
#endif

} // namespace mthesis::nonspectral
