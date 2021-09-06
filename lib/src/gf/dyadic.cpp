#include <mthesis/gf/dyadic.hpp>

namespace mthesis::gf::dyadic
{
namespace free_space
{

DyadC3 G_EJ(const Medium &medium, const VectorR3 &r, const VectorR3 &r_)
{
    using std::complex_literals::operator""i;

    return -1.0i * medium.fd.omega * medium.mu * utils::G_e0(medium, r, r_);

}

DyadC3 G_EM(const Medium &medium, const VectorR3 &r, const VectorR3 &r_)
{
    // Use duality; compare (2.2.35) and (2.2.38) in Jin2015.
    return -utils::G_m0(medium, r, r_);
}

DyadC3 G_HJ(const Medium &medium, const VectorR3 &r, const VectorR3 &r_)
{
    return utils::G_m0(medium, r, r_);
}

DyadC3 G_HM(const Medium &medium, const VectorR3 &r, const VectorR3 &r_)
{
    using std::complex_literals::operator""i;

    // Use duality; compare (2.2.35) and (2.2.38) in Jin2015.
    return -1.0i * medium.fd.omega * medium.eps * utils::G_e0(medium, r, r_);
}

namespace utils {

DyadC3 G_e0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_)
{
    using std::complex_literals::operator""i;

    // See (2.2.36) in Jin2015.

    real k = medium.fd.k_0;
    VectorR3 R_diff = r - r_;
    real R = arma::norm(R_diff);
    cmplx G_0 = std::exp(-1.0i * k * R) / R;
    real x = r(0);
    real y = r(1);
    real z = r(2);
    real x_ = r_(0);
    real y_ = r_(1);
    real z_ = r_(2);

    DyadC3 G;

    G(1, 1) = G_0 * ( -1.0i * std::pow(R, 3) * k - std::pow(R,2) * (std::pow(k, 2) * std::pow(x - x_, 2) + 1) + 3.0i * R * k * std::pow(x - x_, 2) + 3 * std::pow(x - x_, 2) ) / std::pow(R, 4);

    G(2, 1) = G_0 * (x - x_) * (y - y_) * (-std::pow(R, 2) * std::pow(k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(3, 1) = G_0 * (x - x_) * (z - z_) * (-std::pow(R, 2) * std::pow(k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(1, 2) = G_0 * (x - x_) * (y - y_) * (-std::pow(R, 2) * std::pow(k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(2, 2) = G_0 * ( -1.0i * std::pow(R, 3) * k - std::pow(R, 2) * (std::pow(k, 2) * std::pow(y - y_, 2) + 1) + 3.0i * R * k * std::pow(y - y_, 2) + 3 * std::pow(y - y_, 2) ) / std::pow(R, 4);

    G(3, 2) = G_0 * (y - y_) * (z - z_) * (-std::pow(R, 2) * std::pow(k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(1, 3) = G_0 * (x - x_) * (z - z_) * (-std::pow(R, 2) * std::pow(k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(2, 3) = G_0 * (y - y_) * (z - z_) * (-std::pow(R, 2) * std::pow(k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(3, 3) = G_0 * (-1.0i * std::pow(R, 3) * k - std::pow(R, 2) * (std::pow(k, 2) * std::pow(z - z_, 2) + 1) + 3.0i * R * k * std::pow(z - z_, 2) + 3 * std::pow(z - z_, 2) ) / std::pow(R, 4);

    G /= std::pow(k, 2);
    G += arma::diagmat(arma::ones(3)) * G_0;
    G /= 4.0 * M_PI;

    return G;
}

DyadC3 G_m0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_)
{
    using std::complex_literals::operator""i;

    // See (2.2.37) in Jin2015.

    real k = medium.fd.k_0;
    VectorR3 R_diff = r - r_;
    real R = arma::norm(R_diff);
    cmplx G_0 = std::exp(-1.0i * k * R) / R;
    real x = r(0);
    real y = r(1);
    real z = r(2);
    real x_ = r_(0);
    real y_ = r_(1);
    real z_ = r_(2);

    arma::cx_vec3 grad_G;

    grad_G(0) = -G_0 * (x - x_) * (1.0i * R * k + 1.0) / std::pow(R, 2);
    grad_G(1) = -G_0 * (y - y_) * (1.0i * R * k + 1.0) / std::pow(R, 2);
    grad_G(2) = -G_0 * (z - z_) * (1.0i * R * k + 1.0) / std::pow(R, 2);

    DyadC3 G;

    G.col(0) = arma::cross(grad_G, arma::cx_vec3 {1, 0, 0});
    G.col(1) = arma::cross(grad_G, arma::cx_vec3 {1, 0, 0});
    G.col(2) = arma::cross(grad_G, arma::cx_vec3 {1, 0, 0});

    G /= 4.0 * M_PI;

    return G;
}

}

// namespace utils

} // namespace free_space

namespace layered_media
{

namespace utils {

// TODO

} // namespace utils

} // namespace layered_media

} // namespace mthesis::gf::dyadic
