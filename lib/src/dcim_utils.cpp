#include <mthesis/dcim_utils.hpp>
#include <mthesis/gpof.hpp>

#include <cassert>

namespace mthesis::dcim {

std::vector<CmplxImg> get_images(const std::vector<ce_vec> &ce_levels, real rho)
{
    std::vector<CmplxImg> img_all;
    for (const auto &ce_level : ce_levels)
    {
        for (const auto &ce : ce_level)
        {
            auto r = utils::calc_r(rho, ce.exp);
            img_all.emplace_back(ce.amp, r);
        }
    }

    return img_all;
}

cmplx get_spatial_gf(const std::vector<ce_vec> &ce_levels, real rho, real k_0)
{
    using std::complex_literals::operator""i;

    auto img_all = get_images(ce_levels, rho);

    cmplx val = 0;
    for (const auto &img : img_all)
        val += img.amp * std::exp(-1.0i * k_0 * img.r) / img.r;

    val /= 2.0 * M_PI;

    return val;
}

namespace utils {

std::vector<cmplx> get_k_rho_vals(const std::vector<cmplx> &k_z_vals,
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

    for (size_t n = 0; n < k_z_vals.size(); n++)
        k_rho_vals[n] = calc_k_rho(k_z_vals[n]);

    return k_rho_vals;
}

SamplingPath::SamplingPath(real k_0,
                           const std::vector<cmplx> &k_z_vals,
                           real d_t)
    : k_z_vals(k_z_vals),
      k_rho_vals(get_k_rho_vals(k_z_vals, k_0)),
      d_t(d_t),
      N(k_z_vals.size())
{
}

real get_k_0(const LayeredMedium &lm)
{
    auto vacuum = Vacuum(lm.fd);

    real k_min = std::numeric_limits<real>::infinity();

    for (const auto &medium : lm.media)
        if (medium.k.real() < k_min)
            k_min = medium.k.real();

    assert(vacuum.k == k_min);

    return k_min;
}

real find_k_max(const LayeredMedium &lm)
{
    real k_max = -std::numeric_limits<real>::infinity();

    for (const auto &medium : lm.media)
        if (medium.k.real() > k_max)
            k_max = medium.k.real();

    return k_max;
}

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
    if (r.real() < 0.0)
        r *= -1.0;

    return r;
}

cmplx eval_fun(const ce_vec &ce, cmplx k_z)
{
    cmplx val = 0;
    for (const auto &e : ce)
        val += e.amp * std::exp(-e.exp * k_z);

    assert(std::isfinite(val.real()));
    assert(std::isfinite(val.imag()));

    return val;
}

std::vector<cmplx> get_y(const std::function<cmplx (cmplx)> &G,
                         const std::vector<ce_vec> &ce_levels,
                         const std::vector<SamplingPath> &sp, int lev)
{
    using std::complex_literals::operator""i;

    int n_pts = sp[lev].k_rho_vals.size();
    std::vector<cmplx> y(n_pts);

    for (int n = 0; n < n_pts; n++)
    {
        y[n] = G(sp[lev].k_rho_vals[n]) * (1.0i * sp[lev].k_z_vals[n]);

        for (int prev_lev = 0; prev_lev < lev; prev_lev++)
            y[n] -= eval_fun(ce_levels[prev_lev], sp[lev].k_z_vals[n]);

        assert(std::isfinite(y[n].real()));
        assert(std::isfinite(y[n].imag()));
    }

    return y;
}

std::vector<ce_vec> dcim_main_algo(const std::function<cmplx (cmplx)> &G,
                                   const std::vector<SamplingPath> &sp,
                                   const std::vector<ct_fun> &ct_funs)
{
    assert(sp.size() == ct_funs.size());
    int n_lev = sp.size();
    std::vector<ce_vec> ce_levels(n_lev);

    for (int lev = 0; lev < n_lev; lev++)
    {
        auto y = get_y(G, ce_levels, sp, lev);
        auto ce_gpof = gpof::gpof(y, sp[lev].d_t);
        ce_levels[lev] = ct_funs[lev](ce_gpof);
    }
    return ce_levels;
}

} // namespace utils

} // namespace mthesis::dcim
