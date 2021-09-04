#ifndef MTHESIS_SCALAR_GF_HPP
#define MTHESIS_SCALAR_GF_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/tlgf.hpp>
#include <mthesis/sommerfeld_integrals.hpp>

// Generic scalar Green's functions.
namespace mthesis::scalargf {

namespace freespace {

cmplx G_0(const Medium &m, const VectorR3 &r, const VectorR3 &r_)
{
    using std::complex_literals::operator""i;
    real R = arma::norm(r - r_);
    return std::exp(-1.0i * m.k * R) / (4.0 * M_PI * R);
}

} // namespace freespace

namespace layeredmedia {

cmplx generic_spectral(const LayeredMedium &lm,
                          real z,
                          real z_,
                          cmplx k_rho,
                          EmMode mode,
                          bool direct_term)
{
    using std::complex_literals::operator""i;

    auto n = lm.identify_layer(z_);
    bool dual_solution = false;
    const auto d = tlgf::utils::Internals(lm, k_rho, mode, dual_solution);
    auto val = tlgf::utils::V_i_base_generic(lm, z, z_, d, direct_term);

    // Sommerfeld identity + reflected term.
    return val / (2.0i * d.k_z[n]);
}

SommerfeldIntegral get_sommerfeld_integral(const LayeredMedium &lm,
                                           real nu,
                                           EmMode mode,
                                           bool direct_term,
                                           SiParams si_params = SiParams())
{
    auto f = [=](real z, real z_, cmplx k_rho)
    {
        return generic_spectral(lm, z, z_, k_rho, mode, direct_term);
    };

    return SommerfeldIntegral(f, nu, lm, si_params);
}

} // namespace layeredmedia

} // namespace mthesis::scalargf

#endif