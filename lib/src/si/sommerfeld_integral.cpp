#include <mthesis/si/sommerfeld_integral.hpp>
#include <mthesis/si/partition_extrapolation.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <complex_bessel.h>

namespace mthesis {

SommerfeldIntegral::SommerfeldIntegral(SpectralGF f, real nu,
                                       const LayeredMedium &lm)
    : f(f),
      nu(nu),
      lm(lm)
{}

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

} // namespace mthesis
