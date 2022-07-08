#include <mthesis/fmm/gegenbauer.hpp>

#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/legendre.hpp>

namespace mthesis::fmm {

cmplx compute_gegenbauer_pw(const Params &params,
                            const VectorR3 &r,
                            const VectorR3 &r_,
                            unsigned L)
{
    using std::complex_literals::operator""i;

    Group src_group(params, identify_group(params, r_));
    Group obs_group(params, identify_group(params, r));

    // Vector X.
    VectorR3 X = obs_group.r_center - src_group.r_center;
    real X_norm = arma::norm(X);
    VectorR3 X_hat = X / X_norm;
    real kX = params.fd.k_0 * X_norm;

    // Vector d.
    VectorR3 d = r - r_ - X;

    // Precompute plane wave spectrum.
    EwaldSphere es(L);

    // Precompute argument of P_l for all k_hat.
    std::vector<real> cos_theta(es.k_hat.size());
    for (size_t k = 0; k < es.k_hat.size(); k++) {
        cos_theta[k] = arma::dot(es.k_hat[k], X_hat);
    }

    // Evaluate the Gegenbauer sum.
    cmplx val;
    for (unsigned l = 0; l <= L; l++) {
        // Precompute integrand function.
        f_of_k_hat f(es.k_hat.size());
        for (size_t k = 0; k < es.k_hat.size(); k++) {
            f[k] = exp(-1.0i * params.fd.k_0 * arma::dot(es.k_hat[k], d)) *
                    boost::math::legendre_p(l, cos_theta[k]);
        }

        val += pow(-1.0i, l) * (2.0 * l + 1.0) *
                boost::math::sph_hankel_2(l, kX) * es.integrate(f);
    }

    val *= -1.0i * params.fd.k_0 / (4.0 * M_PI);

    val /= 4.0 * M_PI;	// Reference GF has the 4pi in its denominator.

    return val;
}

} // namespace mthesis::fmm
