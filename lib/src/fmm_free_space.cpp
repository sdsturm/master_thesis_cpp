#include <mthesis/fmm_free_space.hpp>

#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/legendre.hpp>

namespace mthesis::fmm::freespace {

unsigned calc_L(const Params &params)
{
    real D = std::sqrt(3.0) * params.w;
    real kD = params.fd.k_0 * D;

    // See (25) in Coifman1993.
    real L = kD + 10.0 * std::log(kD + M_PI);

    return std::ceil(L);
}

f_of_k_hat calc_ff(const FrequencyDomain &fd,
                   const EwaldSphere &es,
                   const Group &g,
                   const VectorR3 &r)
{
    using std::complex_literals::operator""i;
    VectorR3 R = r - g.r_center;
    std::vector<cmplx> V_of_k_hat(es.k_hat.size());
    for (size_t k = 0; k < es.k_hat.size(); k++) {
        V_of_k_hat[k] = std::exp(1.0i * fd.k_0 * arma::dot(es.k_hat[k], R));
    }

    return V_of_k_hat;
}

std::vector<f_of_k_hat> calc_all_ff(const FrequencyDomain &fd,
                                    const EwaldSphere &es,
                                    const std::vector<Group> &groups,
                                    const std::vector<VectorR3> &pts)
{
    size_t K = es.k_hat.size();
    std::vector<f_of_k_hat> ff_all(pts.size(), f_of_k_hat(K));

    for (const auto &g : groups) {
        for (const auto &n : g.pts_idx) {
            ff_all[n] = calc_ff(fd, es, g, pts[n]);
        }
    }

    return ff_all;
}

f_of_k_hat calc_top(unsigned L,
                    const FrequencyDomain &fd,
                    const EwaldSphere &es,
                    const Group &src_group,
                    const Group &obs_group)
{
    using std::complex_literals::operator""i;

    VectorR3 X = obs_group.r_center - src_group.r_center;
    real X_norm = arma::norm(X);
    VectorR3 X_hat = X / X_norm;

    real hankel_arg = fd.k_0 * X_norm;
    std::vector<cmplx> hankel_term(L + 1);
    for (unsigned l = 0; l <= L; l++) {
        hankel_term[l] = std::pow(1.0i, l) * (2.0*l + 1.0) *
                boost::math::sph_hankel_1(l, hankel_arg);
    }

    std::vector<real> cos_theta(es.k_hat.size());
    for (size_t k = 0; k < es.k_hat.size(); k++) {
        cos_theta[k] = arma::dot(es.k_hat[k], X_hat);
    }

    std::vector<cmplx> top(es.k_hat.size());
    cmplx prefactor = 1.0i * fd.k_0 / (4.0 * M_PI);
    for (size_t k = 0; k < es.k_hat.size(); k++) {
        for (unsigned l = 0; l < L; l++) {
            top[k] += hankel_term[l] * boost::math::legendre_p(l, cos_theta[k]);
        }
        top[k] *= prefactor;
    }

    return top;
}

std::vector<f_of_k_hat> calc_all_top(unsigned L,
                                     const FrequencyDomain &fd,
                                     const EwaldSphere &es,
                                     const std::vector<Group> &src_groups,
                                     const std::vector<Group> &obs_groups)
{
    size_t M_ = src_groups.size();
    size_t M = obs_groups.size();
    size_t K = es.k_hat.size();
    using f_of_k_hat = std::vector<cmplx>;
    std::vector<f_of_k_hat> top_all(M_ * M, f_of_k_hat(K));

    size_t i = 0;
    for (const auto &sg : src_groups) {
        for (const auto &og : obs_groups) {
            top_all[i++] = calc_top(L, fd, es, sg, og);
        }
    }

    return top_all;
}

FMM::FMM(const Params &params,
         const std::vector<VectorR3> &src_pts,
         const std::vector<VectorR3> &obs_pts)
    : fd(params.fd),
      w(params.w),
      src_pts(src_pts),
      obs_pts(obs_pts),
      src_groups(build_groups(params, src_pts)),
      obs_groups(build_groups(params, obs_pts)),
      L(calc_L(params)),
      es(EwaldSphere(L)),
      src_ff_all(calc_all_ff(params.fd, es, src_groups, src_pts)),
      obs_ff_all(calc_all_ff(params.fd, es, obs_groups, obs_pts)),
      top_all(calc_all_top(L, params.fd, es, src_groups, obs_groups))
{
    check_group_separation(src_groups, obs_groups, L, fd);
}

std::vector<cmplx> FMM::calc_product(const std::vector<cmplx> &I) const
{
    size_t M_ = src_groups.size();
    size_t M = obs_groups.size();
    size_t K = es.k_hat.size();

    // Step 1: aggregation weighted with source amplitudes.
    std::vector<f_of_k_hat> s(M_, f_of_k_hat(K));
    for (size_t m_ = 0; m_ < M_; m_++) {
        for (const auto &n : src_groups[m_].pts_idx) {
            for (size_t k = 0; k < K; k++) {
                s[m_][k] += std::conj(src_ff_all[n][k]) * I[n];
            }
        }
    }

    // Step 2: translation to observation groups.
    std::vector<f_of_k_hat> g(M, f_of_k_hat(K));
    for (size_t m = 0; m < M; m++) {
        for (size_t m_ = 0; m_ < M_; m_++) {
            size_t top_ind = m_ * M + m; // See loop nesting in calc_all_top.
            for (size_t k = 0; k < K; k++) {
                g[m][k] += top_all[top_ind][k] * s[m_][k];
            }
        }
    }

    // Step 3: disaggregation at observation groups to f(k_hat) and integration.
    std::vector<cmplx> V(obs_pts.size());
    for (size_t m = 0; m < M; m++) {
        for (const auto &n : obs_groups[m].pts_idx) {
            f_of_k_hat f(K);
            for (size_t k = 0; k < K; k++) {
                f[k] = obs_ff_all[n][k] * g[m][k];
            }
            V[n] = es.integrate(f);
        }
    }

    return V;
}

} // namespace mthesis::fmm::freespace
