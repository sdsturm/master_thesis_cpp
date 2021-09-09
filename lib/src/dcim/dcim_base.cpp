#include <mthesis/dcim/dcim_base.hpp>
#include <mthesis/dcim/gpof.hpp>

#include <cassert>

namespace mthesis::dcim {

SamplingPath::SamplingPath(real k_0,
                           const std::vector<cmplx> &k_z_vals,
                           real d_t)
    : k_z_vals(k_z_vals),
      k_rho_vals(get_k_rho_vals(k_z_vals, k_0)),
      d_t(d_t)
{
    assert(d_t > 0.0);
}

std::vector<cmplx> SamplingPath::get_k_rho_vals(
        const std::vector<cmplx> &k_z_vals,
        real k_0)
{
    auto calc_k_rho = [&](cmplx k_z)
    {
        cmplx k_rho = std::sqrt(std::pow(k_0, 2) - std::pow(k_z, 2));

        if (k_rho.imag() < 0.0)
            k_rho *= -1.0;

        return k_rho;
    };

    std::vector<cmplx> k_rho_vals(k_z_vals.size());

    for (size_t n = 0; n < k_z_vals.size(); n++) {
        k_rho_vals[n] = calc_k_rho(k_z_vals[n]);
    }

    return k_rho_vals;
}

CmplxImg::CmplxImg(cmplx amplitude, cmplx r) : amp(amplitude), r(r)
{
    assert(r.real() >= 0);
}

DCIM::DCIM(const SommerfeldIntegral &si) : si(si)
{
    const auto &lm = si.get_lm();

    // Check that the upper layer is semi-infinite and has vacuum parameters.
    assert(std::isinf(lm.d.back()));
    assert(1.0 == lm.media.back().eps_r.real() &&
           0.0 == lm.media.back().eps_r.imag());
    assert(1.0 == lm.media.back().mu_r.real() &&
           0.0 == lm.media.back().mu_r.imag());

    k_0 = lm.fd.k_0;

    k_max = -std::numeric_limits<real>::infinity();
    for (const auto &medium : lm.media) {
        if (medium.k.real() > k_max) {
            k_max = medium.k.real();
        }
    }
}

DCIM::~DCIM() {} // Pure virtual destructor needs to have a definition!

std::vector<CeVec> DCIM::get_exponentials(real z, real z_) const
{
    // Main algorithm.

    auto G = [=](cmplx k_rho) { return si.eval_spectral_gf(z, z_, k_rho); };

    int n_levels = sampling_paths.size();
    std::vector<CeVec> ce_vecs_levels(n_levels);

    for (int lev = 0; lev < n_levels; lev++) {
        auto y = utils::get_y(G, ce_vecs_levels, sampling_paths, lev);
        auto ce_gpof = gpof::gpof(y, sampling_paths[lev].d_t);
        ce_vecs_levels[lev] = ct_funs[lev](ce_gpof);
    }

    return ce_vecs_levels;
}

std::vector<CmplxImg> DCIM::get_images(const std::vector<CeVec> ce_vecs_levels,
                                       real rho) const
{
    std::vector<CmplxImg> imgs_all;
    for (const auto &ce_level : ce_vecs_levels) {
        for (const auto &ce : ce_level) {
            auto r = utils::calc_r(rho, ce.exp);
            imgs_all.emplace_back(ce.amp, r);
        }
    }

    return imgs_all;
}

std::vector<CmplxImg> DCIM::get_images(real z, real z_, real rho) const
{
    const auto ce_vecs_levels = get_exponentials(z, z_);
    return get_images(ce_vecs_levels, rho);
}

cmplx DCIM::get_spatial_gf(const std::vector<CeVec> ce_vecs_levels,
                           real rho) const
{
    using std::complex_literals::operator""i;

    auto img_all = get_images(ce_vecs_levels, rho);

    cmplx val = 0;
    for (const auto &img : img_all) {
        val += img.amp * std::exp(-1.0i * k_0 * img.r) / img.r;
    }

    val /= 2.0 * M_PI;

    return val;
}

cmplx DCIM::get_spatial_gf(real z, real z_, real rho) const
{
    const auto ce_vecs_levels = get_exponentials(z, z_);
    return get_spatial_gf(ce_vecs_levels, rho);
}

cmplx DCIM::get_spatial_gf(const VectorR3 &r, const VectorR3 &r_) const
{
    LayeredMediumCoords coords(r, r_);
    return get_spatial_gf(coords.z, coords.z_, coords.rho);
}



// Hide implementation details in nested namespace utils.
namespace utils {

real get_d_t(real T, int N)
{
    assert(N >= 2);
    return T / (N - 1);
}

cmplx calc_r(real rho, cmplx alpha)
{
    using std::complex_literals::operator""i;

    cmplx arg = std::pow(rho, 2) + std::pow(-1.0i * alpha, 2);
    cmplx r = std::sqrt(arg);

    // See Michalski2007a p. 3 for the branch of the square root.
    if (r.real() < 0.0) {
        r *= -1.0;
    }

    return r;
}

cmplx eval_fun(const CeVec &ce, cmplx k_z)
{
    cmplx val = 0;
    for (const auto &e : ce) {
        val += e.amp * std::exp(-e.exp * k_z);
    }

    assert(std::isfinite(val.real()));
    assert(std::isfinite(val.imag()));

    return val;
}


std::vector<cmplx> get_y(const std::function<cmplx(cmplx)> &G,
                         const std::vector<CeVec> &ce_levels,
                         const std::vector<SamplingPath> &sampling_paths,
                         int level)
{
    using std::complex_literals::operator""i;

    int n_pts = sampling_paths[level].k_rho_vals.size();
    std::vector<cmplx> y(n_pts);

    for (int n = 0; n < n_pts; n++) {
        y[n] = G(sampling_paths[level].k_rho_vals[n]) *
                (1.0i * sampling_paths[level].k_z_vals[n]);

        for (int prev_level = 0; prev_level < level; prev_level++) {
            y[n] -= eval_fun(ce_levels[prev_level],
                             sampling_paths[level].k_z_vals[n]);
        }

        assert(std::isfinite(y[n].real()));
        assert(std::isfinite(y[n].imag()));
    }

    return y;
}

} // namespace utils

} // namespace mthesis::dcim
