#include <mthesis/gf/dyadic.hpp>

namespace mthesis::gf::dyadic
{
namespace free_space
{

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

    G(1, 1) = G_0 * ( -1.0i * pow(R, 3) * k - pow(R,2) * (pow(k, 2) * pow(x - x_, 2) + 1) + 3.0i * R * k * pow(x - x_, 2) + 3 * pow(x - x_, 2) ) / pow(R, 4);

    G(2, 1) = G_0 * (x - x_) * (y - y_) * (-pow(R, 2) * pow(k, 2) + 3.0i * R * k + 3.0) / pow(R, 4);

    G(3, 1) = G_0 * (x - x_) * (z - z_) * (-pow(R, 2) * pow(k, 2) + 3.0i * R * k + 3.0) / pow(R, 4);

    G(1, 2) = G_0 * (x - x_) * (y - y_) * (-pow(R, 2) * pow(k, 2) + 3.0i * R * k + 3.0) / pow(R, 4);

    G(2, 2) = G_0 * ( -1.0i * pow(R, 3) * k - pow(R, 2) * (pow(k, 2) * pow(y - y_, 2) + 1) + 3.0i * R * k * pow(y - y_, 2) + 3 * pow(y - y_, 2) ) / pow(R, 4);

    G(3, 2) = G_0 * (y - y_) * (z - z_) * (-pow(R, 2) * pow(k, 2) + 3.0i * R * k + 3.0) / pow(R, 4);

    G(1, 3) = G_0 * (x - x_) * (z - z_) * (-pow(R, 2) * pow(k, 2) + 3.0i * R * k + 3.0) / pow(R, 4);

    G(2, 3) = G_0 * (y - y_) * (z - z_) * (-pow(R, 2) * pow(k, 2) + 3.0i * R * k + 3.0) / pow(R, 4);

    G(3, 3) = G_0 * (-1.0i * pow(R, 3) * k - pow(R, 2) * (pow(k, 2) * pow(z - z_, 2) + 1) + 3.0i * R * k * pow(z - z_, 2) + 3 * pow(z - z_, 2) ) / pow(R, 4);

    G /= pow(k, 2);
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

    grad_G(0) = -G_0 * (x - x_) * (1.0i * R * k + 1.0) / pow(R, 2);
    grad_G(1) = -G_0 * (y - y_) * (1.0i * R * k + 1.0) / pow(R, 2);
    grad_G(2) = -G_0 * (z - z_) * (1.0i * R * k + 1.0) / pow(R, 2);

    DyadC3 G;

    G.col(0) = arma::cross(grad_G, arma::cx_vec3 {1, 0, 0});
    G.col(1) = arma::cross(grad_G, arma::cx_vec3 {1, 0, 0});
    G.col(2) = arma::cross(grad_G, arma::cx_vec3 {1, 0, 0});

    G /= 4.0 * M_PI;

    return G;
}

} // namespace utils

} // namespace free_space

namespace layered_media
{

} // namespace layered_media

} // namespace mthesis::gf::dyadic
