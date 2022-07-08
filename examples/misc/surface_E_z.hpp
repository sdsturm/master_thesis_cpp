#ifndef MTHESIS_SURFACE_E_Z_HPP
#define MTHESIS_SURFACE_E_Z_HPP

#include <mthesis.hpp>

#include <cassert>

namespace mthesis {

// Get the Sommerfeld integral (46) in Michalski2016b.
inline SommerfeldIntegral get_E_z_surf_VED_si(const LayeredMedium &lm,
                                              real h)
{
    assert(h > 0);
    using std::complex_literals::operator""i;

    auto f = [=](cmplx k_rho, real z, real z_, RiemannSheet sheet)
    {
        tlgf::TLGFParams p(lm, EmMode::TM, true);
        tlgf::utils::Internals d(p, k_rho, false, RiemannSheet::I);

        const Medium &l1 = lm.media.back();
        return -1.0i * l1.eta * l1.k * pow(k_rho, 2) / pow(l1.k, 2) *
                (1.0 - d.Gamma_d.back()) / (2.0i * d.k_z.back()) *
                exp(-1.0i * d.k_z.back() * h);
    };

    return SommerfeldIntegral(f, 0, lm);
}

} // namespace mthesis

#endif // MTHESIS_SURFACE_E_Z_HPP
