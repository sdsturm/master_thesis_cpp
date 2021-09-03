#include <mthesis/fmm.hpp>

#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <random>
#include <cassert>

namespace mthesis::fmm {

Params::Params(const FrequencyDomain &fd, real w) : fd(fd), w(w)
{
    assert(w > 0.0);
}

VectorR3 group_center(const Params &params, const multiindex &mi)
{
    VectorR3 r_center;
    for (unsigned i = 0; i < 3; i++)
        r_center(i) = params.w * (mi(i) + 0.5);

    return r_center;
}

multiindex identify_group(const Params &params, const VectorR3 &r)
{
    multiindex mi;
    for (size_t i = 0; i < 3; i++)
       mi(i) = std::floor((r(i) / params.w));

    return mi;
}

std::vector<VectorR3> rand_pts_in_group(const Params &params,
                                        const multiindex &mi,
                                        unsigned N)
{
    auto r_center = group_center(params, mi);
    std::vector<VectorR3> pts(N, r_center);

    real d = params.w / 2.0;

    std::random_device dev;
    std::mt19937 gen(dev());
    std::uniform_real_distribution<real> dist(-d , d);

    auto rand_rel_point = [&]()
    { return VectorR3 {dist(gen), dist(gen), dist(gen)}; };

    for (auto r : pts)
        r += rand_rel_point();

    return pts;
}

Group::Group(const Params &params, const multiindex &mi)
    : mi(mi), r_center(group_center(params, mi))
{
}

void Group::add_point_index(size_t index)
{
    this->pts_idx.push_back(index);
}

std::vector<Group> build_groups(const Params &params,
                                const std::vector<VectorR3> &pts)
{
    std::vector<Group> glist;

    for (size_t n = 0; n < pts.size(); n++)
    {
        auto mi = identify_group(params, pts[n]);
        bool group_exists = false;
        for (auto g : glist)
        {
            if (arma::all(mi == g.mi))
            {
                group_exists = true;
                g.add_point_index(n);
                break;
            }
        }
        if (!group_exists)
        {
            glist.push_back(Group(params, mi));
            glist.back().add_point_index(n);
        }
    }

    return glist;
}

GroupSeparationError::GroupSeparationError(const Group &sg, const Group &og)
    : sg(sg), og(og)
{
}

void GroupSeparationError::print_message() const
{
    std::cout << "The following groups violate the FMM separation criterion:\n";
    std::cout << "Source group: multiindex = \n" << sg.mi << "\n";
    std::cout << "Observation group: multiindex = \n" << og.mi << "\n";
}

void check_separation(const std::vector<Group> &src_groups,
                      const std::vector<Group> &obs_groups,
                      unsigned L,
                      const FrequencyDomain &fd)
{
    for (const auto &sg : src_groups)
    {
        for (const auto &og : obs_groups)
        {
            real R = arma::norm(og.r_center - sg.r_center);
            real factor = 1.0;	// TODO: play around with this.
            auto L_max = factor * fd.lambda_0 * R;
            if (L > L_max)
                throw GroupSeparationError(sg, og);
        }
    }
}

EwaldSpere::EwaldSpere(unsigned L)
    : cos_theta_weights(L)
{
    // Get nodes and weights for Gauss-Legendre quadrature in [-1, 1].
    std::vector<real> cos_thets_nodes(cos_theta_weights.size());
    auto t = gsl_integration_glfixed_table_alloc(cos_theta_weights.size());
    for (size_t n = 0; n < cos_theta_weights.size(); n++)
    {
        gsl_integration_glfixed_point(-1.0, 1.0, n,
                                      cos_thets_nodes.data() + n,
                                      cos_theta_weights.data() + n,
                                      t);
    }
    gsl_integration_glfixed_table_free(t);

    // Get nodes in phi for [0, 2pi).
    size_t n_phi_nodes = 2 * L;
    real phi_step = 2.0 * M_PI / n_phi_nodes;
    std::vector<real> phi_nodes(n_phi_nodes);
    for (size_t n = 0; n < n_phi_nodes; n++)
    {
        phi_nodes[n] = n * phi_step;
    }

    this->prefactor = phi_step;

    k_hat.resize(cos_theta_weights.size() * n_phi_nodes);
    size_t i = 0;
    for (const auto &phi : phi_nodes)
    {
        for (const auto &cos_theta : cos_thets_nodes)
        {
            real sin_theta = std::sqrt(1.0 - std::pow(cos_theta, 2.0));
            real k_x = std::cos(phi) * sin_theta;
            real k_y = std::sin(phi) * sin_theta;
            real k_z = cos_theta;
            k_hat[i] = VectorR3{k_x, k_y, k_z};
            i++;
        }
    }
}

cmplx EwaldSpere::integrate(const std::vector<cmplx> &f_of_k_hat) const
{
    assert(k_hat.size() == f_of_k_hat.size());

    cmplx val = 0;
    for (size_t n = 0; n < k_hat.size(); n++)
        val += cos_theta_weights[n] * f_of_k_hat[n];

    return prefactor * val;
}

unsigned calc_L(const Params &params)
{
    real D = std::sqrt(3.0) * params.w;
    real kD = params.fd.k_0 * D;

    // See (25) in Coifman1993.
    real L = kD + 10.0 * std::log(kD + M_PI);

    return std::ceil(L);
}

std::vector<cmplx> calc_ff(const FrequencyDomain &fd,
                           const EwaldSpere &es,
                           const Group &g,
                           const VectorR3 &r)
{
    using std::complex_literals::operator""i;
    VectorR3 R = r - g.r_center;
    std::vector<cmplx> V_of_k_hat(es.k_hat.size());
    for (size_t k = 0; k < es.k_hat.size(); k++)
        V_of_k_hat[k] = std::exp(1.0i * fd.k_0 * arma::dot(es.k_hat[k], R));

    return V_of_k_hat;
}

std::vector<cmplx> calc_top(unsigned L,
                            const FrequencyDomain &fd,
                            const EwaldSpere &es,
                            const Group &src_group,
                            const Group &obs_group)
{
    using std::complex_literals::operator""i;

    VectorR3 X = obs_group.r_center - src_group.r_center;
    real X_norm = arma::norm(X);
    VectorR3 X_hat = X / X_norm;

    real hankel_arg = fd.k_0 * X_norm;
    std::vector<cmplx> hankel_term(L + 1);
    for (unsigned l = 0; l <= L; l++)
    {
        hankel_term[l] = std::pow(1.0i, l) * (2.0*l + 1.0) *
                boost::math::sph_hankel_1(l, hankel_arg);
    }

    std::vector<real> cos_theta(es.k_hat.size());
    for (size_t k = 0; k < es.k_hat.size(); k++)
    {
        cos_theta[k] = arma::dot(es.k_hat[k], X_hat);
    }

    std::vector<cmplx> top(es.k_hat.size());
    for (size_t k = 0; k < es.k_hat.size(); k++)
    {
        for (unsigned l = 0; l < L; l++)
        {
            top[k] += hankel_term[l] * boost::math::legendre_p(l, cos_theta[k]);
        }
        top[k] *= 1.0i * fd.k_0 / (4.0 * M_PI);
    }

    return top;
}

FreeSpaceFMM::FreeSpaceFMM(const Params &params,
                           const std::vector<VectorR3> &src_pts,
                           const std::vector<VectorR3> &obs_pts)
    : fd(params.fd),
      w(params.w),
      src_pts(src_pts),
      obs_pts(obs_pts),
      src_groups(build_groups(params, src_pts)),
      obs_groups(build_groups(params, obs_pts)),
      L(calc_L(params)),
      es(EwaldSpere(L))
{
}

} // namespace mthesis::fmm
