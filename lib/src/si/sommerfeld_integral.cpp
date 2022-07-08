#include <mthesis/si/sommerfeld_integral.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <complex_bessel.h>

namespace mthesis {

SommerfeldIntegral::SommerfeldIntegral(SpectralGF f, real nu,
                                       const LayeredMedium &lm)
    : f(f),
      nu(nu),
      lm(lm)
{}

cmplx SommerfeldIntegral::eval_spectral_gf(cmplx k_rho,
                                           real z,
                                           real z_,
                                           RiemannSheet sheet) const
{
    return f(k_rho, z, z_, sheet);
}

cmplx SommerfeldIntegral::eval_integrand_sip(real k_rho,
                                             real rho,
                                             real z,
                                             real z_,
                                             RiemannSheet sheet) const
{
    return eval_spectral_gf(k_rho, z, z_, sheet) * k_rho *
            boost::math::cyl_bessel_j(nu, rho * k_rho);
}

cmplx SommerfeldIntegral::eval_integrand_sip(cmplx k_rho,
                                             real rho,
                                             real z,
                                             real z_,
                                             RiemannSheet sheet) const
{
    return eval_spectral_gf(k_rho, z, z_, sheet) * k_rho *
            sp_bessel::besselJ(nu, rho * k_rho);

}

} // namespace mthesis
