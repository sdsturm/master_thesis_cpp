#include <mthesis/sommerfeld_integrals.hpp>
#include <mthesis/partition_extrapolation.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <complex_bessel.h>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cmath>
#include <cassert>

namespace mthesis {

SiParams::SiParams()
    : alpha(std::numeric_limits<real>::quiet_NaN()),
      zeta(std::numeric_limits<real>::quiet_NaN()),
      identify_singularities(false)
{
    assert( (std::isnan(alpha) && std::isnan(zeta)) ||
            (std::isfinite(alpha) && std::isfinite(zeta)) );

    if (std::isfinite(zeta))
        assert(zeta >= 0.0);

}

void SiParams::set_alpha(real alpha)
{
    assert(std::isfinite(alpha));
    this->alpha = alpha;
}

void SiParams::set_zeta(real zeta)
{
    assert(zeta >= 0.0);
    this->zeta = zeta;
}

void SiParams::set_identify_singularities(bool identify_singularities)
{
    this->identify_singularities = identify_singularities;

}

SommerfeldIntegral::SommerfeldIntegral(spectral_gf f,
                                       real nu,
                                       const LayeredMedium &lm,
                                       SiParams params)
    : f(f),
      nu(nu),
      lm(lm),
      params(params),
      bp(),		// not determined by default
      swp()		// not determined by default
{
    assert(nu >= 0.0);

    if (params.identify_singularities)
    {
        this->bp = get_branch_points(lm);

        assert(false);	// Not implemented yet.
        this->swp = identify_swp();
    }
}

cmplx SommerfeldIntegral::eval_spectral_gf(real z, real z_, cmplx k_rho) const
{
    return f(z, z_, k_rho);
}

cmplx SommerfeldIntegral::eval_integrand_sip(real rho, real z, real z_,
                                             real k_rho) const
{
    return f(z, z_, k_rho) * k_rho * boost::math::cyl_bessel_j(nu, rho * k_rho);
}

cmplx SommerfeldIntegral::eval_integrand_sip(real rho, real z, real z_,
                                             cmplx k_rho) const
{
    return f(z, z_, k_rho) * k_rho * sp_bessel::besselJ(nu, rho * k_rho);
}

cmplx SommerfeldIntegral::eval_si_along_sip(real rho, real z, real z_,
                                            pe::Params pe_params) const
{
    if (rho == 0.0 && nu != 0.0)
        return 0.0;

    auto a = get_a();

    auto head = eval_head_ellipsis(rho, z, z_, a);
    auto tail = eval_tail_on_sip(rho, z, z_, a, pe_params);

    return (head + tail) / (2.0 * M_PI);
}

cmplx SommerfeldIntegral::eval_si_along_sip(const VectorR3 &r,
                                            const VectorR3 &r_,
                                            pe::Params pe_params) const
{
    auto coords = LayeredMediumCoords(r, r_);

    return eval_si_along_sip(coords.rho, coords.z, coords.z_, pe_params);
}

cmplx SommerfeldIntegral::eval_integrand_eip(real rho, real z, real z_,
                                             cmplx k_rho) const
{
    return f(z, z_, k_rho) * k_rho * sp_bessel::hankelH2(nu, rho * k_rho);
}

const LayeredMedium &SommerfeldIntegral::get_lm() const
{
    return lm;
}

real SommerfeldIntegral::get_a() const
{
    std::vector<real> k_re(lm.media.size());
    for (size_t n = 0; n < lm.media.size(); n++)
        k_re[n] = lm.media[n].k.real();

    real k_re_max = *std::max_element(k_re.begin(), k_re.end());

    real a = 1.2 * k_re_max;

    return a;
}

real SommerfeldIntegral::calc_indention(real rho, real a) const
{
    real w;
    if (rho > 0.0)
        w = 1.0 / rho;
    else
        w = a / 2.0;

    return w;
}

cmplx SommerfeldIntegral::eval_head_ellipsis(real rho,
                                             real z,
                                             real z_,
                                             real a) const
{
    using std::complex_literals::operator""i;

    real r = a / 2.0;
    real w = calc_indention(rho, a);

    auto f_along_path = [&](real t)
    {
        // For 0 <= t <= 1.
        real arg = M_PI * (1.0 - t);
        auto sin_val = std::sin(arg);
        auto cos_val = std::cos(arg);
        cmplx gamma = r * (1.0 + cos_val) + 1.0i * w * sin_val;
        cmplx gamma_ = r * M_PI * sin_val - 1.0i * w * M_PI * cos_val;
        return eval_integrand_sip(rho, z, z_, gamma) * gamma_;
    };

    constexpr boost::math::quadrature::gauss_kronrod<real, 31> quad;
    auto val = quad.integrate(f_along_path, 0.0, 1.0);

    return val;
}

cmplx SommerfeldIntegral::eval_tail_on_sip(real rho,
                                           real z,
                                           real z_,
                                           real a,
                                           const pe::Params &pe_params) const
{
    auto integrand = [=](real k_rho)
    { return eval_integrand_sip(rho, z, z_, k_rho); };

    cmplx val;
    if (rho > 0.0)
    {
        if (std::isfinite(params.alpha) && std::isfinite(params.zeta))
        {
            return pe::mosig_michalski(integrand,
                                       params.alpha, params.zeta,
                                       nu, rho, a,
                                       pe_params);
        }
        else
        {
            return pe::levin_sidi(integrand, nu, rho, a, pe_params);
        }
    }
    else if (rho == 0 && std::abs(z - z_) > 0.0)
    {
        constexpr real b = std::numeric_limits<real>::infinity();
        constexpr boost::math::quadrature::gauss_kronrod<real, 15> quad;
        val = quad.integrate(integrand, a, b);
    }
    else
    {
        std::cerr << "SI with rho = |z - z_| = 0. Returning NaN.";
        val = std::numeric_limits<real>::quiet_NaN();
    }
    return val;

}

std::vector<cmplx> get_branch_points(const LayeredMedium &lm)
{
    std::vector<cmplx> bp_locations;

    if (std::isinf(lm.d.front()))
        bp_locations.push_back(lm.media.front().k);

    if (std::isinf(lm.d.back()))
        bp_locations.push_back(lm.media.back().k);

    return bp_locations;
}

std::vector<cmplx> identify_swp(/* TODO */)
{
    std::vector<cmplx> sip_locations;

    // TODO

    return sip_locations;
}

} // namespace mthesis
