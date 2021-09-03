#include <mthesis/sgf.hpp>
#include <mthesis/sommerfeld_integrals.hpp>

namespace mthesis::sgf {

cmplx free_space(const Medium &m, const VectorR3 &r, const VectorR3 &r_)
{
    using std::complex_literals::operator""i;
    real R = arma::norm(r - r_);
    return std::exp(-1.0i * m.k * R) / (4.0 * M_PI * R);
}

cmplx lm_generic_spectral(const LayeredMedium &lm,
                          real z,
                          real z_,
                          cmplx k_rho,
                          EmMode mode,
                          bool direct_term)
{
    using std::complex_literals::operator""i;

    auto n = lm.identify_layer(z_);
    bool dual_solution = false;
    const auto d = tlgf::Internals(lm, k_rho, mode, dual_solution);
    auto val = tlgf::V_i_base_generic(lm, z, z_, d, direct_term);

    // Sommerfeld identity + reflected term.
    return val / (2.0i * d.k_z[n]);
}

const si::SpectralGF lm_get_generic_spectral_gf(const LayeredMedium &lm,
                                                const VectorR3 &r,
                                                const VectorR3 &r_,
                                                EmMode mode,
                                                bool direct_term)
{
    auto coords = LayeredMediumCoords(r, r_);

    auto f = [=](cmplx k_rho)
    {
        return lm_generic_spectral(lm, coords.z, coords.z_, k_rho, mode, direct_term);
    };

    const auto gf = si::SpectralGF(f, coords.z, coords.z_, lm, false);

    return gf;
}

cmplx lm_generic_spatial(const LayeredMedium &lm,
                         const VectorR3 &r,
                         const VectorR3 &r_,
                         EmMode mode,
                         bool direct_term)
{
    auto coords = LayeredMediumCoords(r, r_);
    auto gf = lm_get_generic_spectral_gf(lm, r, r_, mode, direct_term);

    real nu = 0;
    si::pe::Params pe_params;
    return si::eval_along_sip(gf, nu, coords.rho, pe_params);
}

} // namespace mthesis::sgf
