#include <mthesis/gf/scalar.hpp>

#include <mthesis/tlgf.hpp>
#include <mthesis/sommerfeld_integrals.hpp>

namespace mthesis::gf::scalar {

namespace free_space {

cmplx G_0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_)
{
    using std::complex_literals::operator""i;
    real R = arma::norm(r - r_);
    return std::exp(-1.0i * medium.k * R) / (4.0 * M_PI * R);
}

cmplx G_0(const Medium &medium, const VectorR3 &r, const VectorC3 &r_)
{
    using std::complex_literals::operator""i;
    cmplx R = cmplx_length(r - r_);
    return std::exp(-1.0i * medium.k * R) / (4.0 * M_PI * R);
}

} // namespace free_space

namespace layered_media {

cmplx generic_spectral_gf(const LayeredMedium &lm,
                          cmplx k_rho,
                          real z,
                          real z_,
                          EmMode mode,
                          bool direct_term)
{
    using std::complex_literals::operator""i;

    cmplx val;
    size_t n = lm.identify_layer(z_);

    // Note: this function is not reciprocal if z and z_ lie in layers with
    // 	     different material paramters.
    //       This is due to the division by possibliy different eps or mu
    //       in both argument cases (z, z_) and (z_, z).
    //       For the case m == n (z and z_ in the same layer) this function
    //       fulfills reciprocity.

    switch (mode) {
    case (EmMode::TM):
        val = tlgf::I_v(lm, k_rho, z, z_, mode, direct_term);
        val /= 1.0i * lm.fd.omega * lm.media[n].eps;
        break;
    case (EmMode::TE):
        val = tlgf::V_i(lm, k_rho, z, z_, mode, direct_term);
        val /= 1.0i * lm.fd.omega * lm.media[n].mu;
        break;
    }

    return val;
}

SommerfeldIntegral get_sommerfeld_integral(const LayeredMedium &lm,
                                           real nu,
                                           EmMode mode,
                                           bool direct_term,
                                           SiParams si_params)
{
    auto f = [=](real z, real z_, cmplx k_rho)
    {
        return generic_spectral_gf(lm, k_rho, z, z_, mode, direct_term);
    };

    return SommerfeldIntegral(f, nu, lm, si_params);
}

} // namespace layered_media

} // namespace mthesis::gf::scalar
