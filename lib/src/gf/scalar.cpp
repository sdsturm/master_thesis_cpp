#include <mthesis/gf/scalar.hpp>

#include <mthesis/tlgf.hpp>

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
                                           bool direct_term,
                                           RiemannSheet s)
{
    auto f = [=](real z, real z_, cmplx k_rho)
    {
        tlgf::TLGFParams p(lm, mode, direct_term, s);
        return tlgf::generic_sgf(k_rho, z, z_, p);
    };

    return SommerfeldIntegral(f, nu, lm);
}

} // namespace layered_media

} // namespace mthesis::gf::scalar
