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

SommerfeldIntegral get_sommerfeld_integral(const LayeredMedium &lm,
                                           real nu,
                                           EmMode mode,
                                           bool direct_term)
{
    auto f = [=](real z, real z_, cmplx k_rho)
    {
        return tlgf::generic_sgf(lm, k_rho, z, z_, mode, direct_term);
    };

    return SommerfeldIntegral(f, nu, lm);
}

} // namespace layered_media

} // namespace mthesis::gf::scalar
