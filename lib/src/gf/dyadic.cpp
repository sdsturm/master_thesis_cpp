#include <mthesis/gf/dyadic.hpp>
#include <mthesis/tlgf.hpp>
#include <mthesis/si/sommerfeld_integral.hpp>
#include <mthesis/si/axial_transmission.hpp>

#include <gsl/gsl_const_mksa.h>

#include <cassert>

namespace mthesis::gf::dyadic {
namespace free_space {

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

// Hide details.
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

    G(0, 0) = G_0 * ( -1.0i * std::pow(R, 3) * k - std::pow(R, 2) * ( std::pow(k * (x - x_), 2) + 1.0 ) + 3.0i * R * k * std::pow(x - x_, 2) + 3.0 * std::pow(x - x_, 2) ) / std::pow(R, 4);

    G(1, 0) = G_0 * (x - x_) * (y - y_) * (-std::pow(R * k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(2, 0) = G_0 * (x - x_) * (z - z_) * (-std::pow(R * k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(0, 1) = G_0 * (x - x_) * (y - y_) * (-std::pow(R * k, 2) + 3.0i * R * k + 3.0) / std::pow(R, 4);

    G(1, 1) = G_0 * ( -1.0i * std::pow(R, 3) * k - std::pow(R, 2) * ( std::pow(k * (y - y_), 2) + 1.0 ) + 3.0i * R * k * std::pow(y - y_, 2) + 3.0 * std::pow(y - y_, 2) ) / std::pow(R, 4);

    G(2, 1) = G_0 * (y - y_) * (z - z_) * ( -std::pow(R * k, 2) + 3.0i * R * k + 3.0 ) / std::pow(R, 4);

    G(0, 2) = G_0 * (x - x_) * (z - z_) * ( -std::pow(R * k, 2) + 3.0i * R * k + 3.0 ) / std::pow(R, 4);

    G(1, 2) = G_0 * (y - y_) * (z - z_) * ( -std::pow(R * k, 2) + 3.0i * R * k + 3.0 ) / std::pow(R, 4);

    G(2, 2) = G_0 * ( -1.0i * std::pow(R, 3) * k - std::pow(R, 2) * ( std::pow(k * (z - z_), 2) + 1.0 ) + 3.0i * R * k * std::pow(z - z_, 2) + 3.0 * std::pow(z - z_, 2) ) / std::pow(R, 4);

    G /= std::pow(k, 2);

    arma::cx_mat33 I(arma::fill::eye);
    G += I * G_0;

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

    cmplx term = -G_0 * (1.0i * R * k + 1.0) / std::pow(R, 2);
    grad_G(0) = term * (x - x_);
    grad_G(1) = term * (y - y_);
    grad_G(2) = term * (z - z_);

    DyadC3 G;

    G(0, 0) = 0;
    G(1, 0) = grad_G(2);
    G(2, 0) = -grad_G(1);

    G(0, 1) = -grad_G(2);
    G(1, 1) = 0;
    G(2, 1) = grad_G(0);

    G(0, 2) = grad_G(1);
    G(1, 2) = -grad_G(0);
    G(2, 2) = 0;

    G /= 4.0 * M_PI;

    return G;
}

}

// namespace utils

} // namespace free_space

namespace layered_media {

DyadC3 G_EJ(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_)
{
    LMCoords c(r, r_);

    auto si_vals = utils::calc_G_EJ_si_vals(lm, c);
    auto f_ = utils::calc_electric_factor(lm, c.z_);
    auto f = utils::calc_electric_factor(lm, c.z);

    return utils::G_EJ_HM_common(si_vals, c, f, f_);
}

DyadC3 G_HM(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_)
{
    LMCoords c(r, r_);

    auto si_vals = utils::calc_G_HM_si_vals(lm, c);
    auto f_ = utils::calc_magnetic_factor(lm, c.z_);
    auto f = utils::calc_magnetic_factor(lm, c.z);

    return utils::G_EJ_HM_common(si_vals, c, f, f_);
}


DyadC3 G_HJ(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_)
{
    LMCoords c(r, r_);

    auto si_vals = utils::calc_G_HJ_si_vals(lm, c);

    auto factor_e = utils::calc_electric_factor(lm, c.z_);
    auto factor_h = utils::calc_magnetic_factor(lm, c.z);

    DyadC3 G;

    G(0, 0) = utils::G_HJ_xx(si_vals, c);
    G(1, 0) = utils::G_HJ_yx(si_vals, c);
    G(2, 0) = utils::G_HJ_zx(si_vals, c, factor_h);

    G(0, 1) = utils::G_HJ_xy(si_vals, c);
    G(1, 1) = -G(0, 0);
    G(2, 1) = utils::G_HJ_zy(si_vals, c, factor_h);

    G(0, 2) = utils::G_HJ_xz(si_vals, c, factor_e);
    G(1, 2) = utils::G_HJ_yz(si_vals, c, factor_e);
    G(2, 2) = 0;

    return G;
}

DyadC3 G_EM(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_)
{
    // Use reciprocity relation (10) in Michalski2005.
    return -arma::strans(G_HJ(lm, r_, r));
}

// Hide details.
namespace utils {

real get_eta_0()
{
    return sqrt(GSL_CONST_MKSA_VACUUM_PERMEABILITY /
                GSL_CONST_MKSA_VACUUM_PERMITTIVITY);
}

cmplx calc_electric_factor(const LayeredMedium &lm, real z)
{
    using std::complex_literals::operator""i;

    size_t m = lm.identify_layer(z);

    return get_eta_0() / (1.0i * lm.fd.k_0 * lm.media[m].eps_r);
}

cmplx calc_magnetic_factor(const LayeredMedium &lm, real z)
{
    using std::complex_literals::operator""i;

    size_t m = lm.identify_layer(z);

    return 1.0 / (1.0i * lm.fd.k_0 * get_eta_0() * lm.media[m].mu_r);
}

std::vector<cmplx> calc_G_EJ_si_vals(const LayeredMedium &lm, const LMCoords &c)
{
    using namespace mthesis::tlgf;

    auto f0 = [=](real z, real z_, cmplx k_rho)
    {
        return V_i(k_rho, z, z_, TLGFParams(lm, EmMode::TM));
    };

    auto f1 = [=](real z, real z_, cmplx k_rho)
    {
        return V_i(k_rho, z, z_, TLGFParams(lm, EmMode::TE));
    };

    auto f2 = [=](real z, real z_, cmplx k_rho)
    {
        return (f0(z, z_, k_rho) - f1(z, z_, k_rho)) / k_rho;
    };

    auto f3 = [=](real z, real z_, cmplx k_rho)
    {
        return k_rho * V_v(k_rho, z, z_, TLGFParams(lm, EmMode::TM));
    };

    auto f4 = [=](real z, real z_, cmplx k_rho)
    {
        return k_rho * I_i(k_rho, z, z_, TLGFParams(lm, EmMode::TM));
    };

    auto f5 = [=](real z, real z_, cmplx k_rho)
    {
        return pow(k_rho, 2) * I_v(k_rho, z, z_, TLGFParams(lm, EmMode::TM));
    };

    std::vector<SommerfeldIntegral> si;
    si.push_back(SommerfeldIntegral(f0, 0, lm));
    si.push_back(SommerfeldIntegral(f1, 0, lm));
    si.push_back(SommerfeldIntegral(f2, 1, lm));
    si.push_back(SommerfeldIntegral(f3, 1, lm));
    si.push_back(SommerfeldIntegral(f4, 1, lm));
    si.push_back(SommerfeldIntegral(f5, 0, lm));

    using mthesis::si::axial_transmission::eval_si_along_sip;
    std::vector<cmplx> si_vals(6);
    for (size_t n = 0; n < si.size(); n++) {
        si_vals[n] = eval_si_along_sip(si[n], c.rho, c.z, c.z_);
    }

    return si_vals;
}

std::vector<cmplx> calc_G_HM_si_vals(const LayeredMedium &lm, const LMCoords &c)
{
    using namespace mthesis::tlgf;

    auto f0 = [=](real z, real z_, cmplx k_rho)
    {
        return I_v(k_rho, z, z_, TLGFParams(lm, EmMode::TE));
    };

    auto f1 = [=](real z, real z_, cmplx k_rho)
    {
        return I_v(k_rho, z, z_, TLGFParams(lm, EmMode::TM));
    };

    auto f2 = [=](real z, real z_, cmplx k_rho)
    {
        return (f0(z, z_, k_rho) - f1(z, z_, k_rho)) / k_rho;
    };

    auto f3 = [=](real z, real z_, cmplx k_rho)
    {
        return k_rho * I_i(k_rho, z, z_, TLGFParams(lm, EmMode::TE));
    };

    auto f4 = [=](real z, real z_, cmplx k_rho)
    {
        return k_rho * V_v(k_rho, z, z_, TLGFParams(lm, EmMode::TE));
    };

    auto f5 = [=](real z, real z_, cmplx k_rho)
    {
        return pow(k_rho, 2) * V_i(k_rho, z, z_, TLGFParams(lm, EmMode::TE));
    };

    std::vector<SommerfeldIntegral> si;
    si.push_back(SommerfeldIntegral(f0, 0, lm));
    si.push_back(SommerfeldIntegral(f1, 0, lm));
    si.push_back(SommerfeldIntegral(f2, 1, lm));
    si.push_back(SommerfeldIntegral(f3, 1, lm));
    si.push_back(SommerfeldIntegral(f4, 1, lm));
    si.push_back(SommerfeldIntegral(f5, 0, lm));

    using mthesis::si::axial_transmission::eval_si_along_sip;
    std::vector<cmplx> si_vals(6);
    for (size_t n = 0; n < si.size(); n++) {
        si_vals[n] = eval_si_along_sip(si[n], c.rho, c.z, c.z_);
    }

    return si_vals;
}

cmplx G_EJ_HM_xx(const std::vector<cmplx> &si_vals, const LMCoords &c)
{
    cmplx val = -pow(cos(c.phi), 2) * si_vals[0] -
            pow(sin(c.phi), 2) * si_vals[1];

    if (c.rho > 0.0) {
        val += cos(2.0 * c.phi) / c.rho * si_vals[2];
    }
    return val;
}

cmplx G_EJ_HM_xz(const std::vector<cmplx> &si_vals, const LMCoords &c,
                 cmplx f_)
{
    return f_ * cos(c.phi) * si_vals[3];
}

cmplx G_EJ_HM_yx(const std::vector<cmplx> &si_vals, const LMCoords &c)
{
    cmplx val = -(si_vals[0] - si_vals[1]) / 2.0;

    if (c.rho > 0.0) {
        val += si_vals[2] / c.rho;
    }

    return sin(2.0 * c.phi) * val;
}

cmplx G_EJ_HM_yy(const std::vector<cmplx> &si_vals, const LMCoords &c)
{
    cmplx val = -pow(sin(c.phi), 2) * si_vals[0] -
            pow(cos(c.phi), 2) * si_vals[1];

    if (c.rho > 0) {
        val -= cos(2.0 * c.phi) / c.rho * si_vals[2];
    }

    return val;
}

cmplx G_EJ_HM_yz(const std::vector<cmplx> &si_vals, const LMCoords &c,
                 cmplx f_)
{
    return f_ * sin(c.phi) * si_vals[3];
}

cmplx G_EJ_HM_zx(const std::vector<cmplx> &si_vals, const LMCoords &c,
                 cmplx f)
{
    return f * cos(c.phi) * si_vals[4];
}

cmplx G_EJ_HM_zy(const std::vector<cmplx> &si_vals, const LMCoords &c,
                 cmplx f)
{
    return f * sin(c.phi) * si_vals[4];
}

cmplx G_EJ_HM_zz(const std::vector<cmplx> &si_vals, const LMCoords &c,
                 cmplx f, cmplx f_)
{
    cmplx val = f_ * f * si_vals[5];

    if (c.rho == 0.0 && (c.z - c.z_) == 0.0) {
        // A delta-function outside integral seems quite strange...
        val -= f * std::numeric_limits<real>::infinity();
    }

    return val;
}

DyadC3 G_EJ_HM_common(const std::vector<cmplx> &si_vals, const LMCoords &c,
                      cmplx f, cmplx f_)
{
    DyadC3 G;

    G(0, 0) = G_EJ_HM_xx(si_vals, c);
    G(1, 0) = G_EJ_HM_yx(si_vals, c);
    G(2, 0) = G_EJ_HM_zx(si_vals, c, f);

    G(0, 1) = G(1, 0);
    G(1, 1) = G_EJ_HM_yy(si_vals, c);
    G(2, 1) = G_EJ_HM_zy(si_vals, c, f);

    G(0, 2) = G_EJ_HM_xz(si_vals, c, f_);
    G(1, 2) = G_EJ_HM_yz(si_vals, c, f_);
    G(2, 2) = G_EJ_HM_zz(si_vals, c, f, f_);

    return G;
}

std::vector<cmplx> calc_G_HJ_si_vals(const LayeredMedium &lm, const LMCoords &c)
{
    using namespace mthesis::tlgf;

    auto f0 = [=](real z, real z_, cmplx k_rho)
    {
        return I_i(k_rho, z, z_, TLGFParams(lm, EmMode::TE));
    };

    auto f1 = [=](real z, real z_, cmplx k_rho)
    {
        return I_i(k_rho, z, z_, TLGFParams(lm, EmMode::TM));
    };

    auto f2 = [=](real z, real z_, cmplx k_rho)
    {
        return (f0(z, z_, k_rho) - f1(z, z_, k_rho)) / k_rho;
    };

    auto f3 = [=](real z, real z_, cmplx k_rho)
    {
        return k_rho * I_v(k_rho, z, z_, TLGFParams(lm, EmMode::TM));
    };

    auto f4 = [=](real z, real z_, cmplx k_rho)
    {
        return k_rho * V_i(k_rho, z, z_, TLGFParams(lm, EmMode::TE));
    };

    std::vector<SommerfeldIntegral> si;
    si.push_back(SommerfeldIntegral(f0, 0, lm));
    si.push_back(SommerfeldIntegral(f1, 0, lm));
    si.push_back(SommerfeldIntegral(f2, 1, lm));
    si.push_back(SommerfeldIntegral(f3, 1, lm));
    si.push_back(SommerfeldIntegral(f4, 1, lm));

    using mthesis::si::axial_transmission::eval_si_along_sip;
    std::vector<cmplx> si_vals(5);
    for (size_t n = 0; n < si.size(); n++) {
        si_vals[n] = eval_si_along_sip(si[n], c.rho, c.z, c.z_);
    }

    return si_vals;
}

cmplx G_HJ_xx(const std::vector<cmplx> &si_vals, const LMCoords &c)
{
    cmplx val = -(si_vals[0] - si_vals[1]) / 2.0;

    if (c.rho > 0.0) {
        val += si_vals[2] / c.rho;
    }

    return sin(2.0 * c.phi) * val;
}

cmplx G_HJ_xy(const std::vector<cmplx> &si_vals, const LMCoords &c)
{
    cmplx val = pow(cos(c.phi), 2) * si_vals[0] +
            pow(sin(c.phi), 2) * si_vals[1];

    if (c.rho > 0.0) {
        val -= cos(2.0 * c.phi) / c.rho * si_vals[2];
    }

    return val;
}

cmplx G_HJ_xz(const std::vector<cmplx> &si_vals, const LMCoords &c,
              cmplx factor_e_)
{
    return -factor_e_ * sin(c.phi) * si_vals[3];
}

cmplx G_HJ_yx(const std::vector<cmplx> &si_vals, const LMCoords &c)
{
    cmplx val = -pow(sin(c.phi), 2) * si_vals[0] -
            pow(cos(c.phi), 2) * si_vals[1];

    if (c.rho > 0.0) {
        val -= cos(2.0 * c.phi) / c.rho * si_vals[2];
    }

    return val;
}

cmplx G_HJ_yz(const std::vector<cmplx> &si_vals, const LMCoords &c,
              cmplx factor_e_)
{
    return factor_e_ * cos(c.phi) * si_vals[3];
}

cmplx G_HJ_zx(const std::vector<cmplx> &si_vals, const LMCoords &c,
              cmplx factor_h)
{
    return factor_h * sin(c.phi) * si_vals[4];
}

cmplx G_HJ_zy(const std::vector<cmplx> &si_vals, const LMCoords &c,
              cmplx factor_h)
{
    return -factor_h * cos(c.phi) * si_vals[4];
}

} // namespace utils

} // namespace layered_media

} // namespace mthesis::gf::dyadic
