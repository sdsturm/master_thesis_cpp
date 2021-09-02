#include <mthesis/detail/dcim_utils.hpp>
#include <mthesis/detail/gpof.hpp>

#include <cassert>

namespace mthesis::dcim::utils
{
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

        for (const auto e : ce)
            val += e.front() * std::exp(-e.back() * k_z);

        return val;
    }

    std::vector<ce_vec> algo(const si::SpectralGF &gf,
                             const std::vector<SamplingPath> &sp,
                             const std::vector<ct_fun> &ct_funs)
    {
        using std::complex_literals::operator""i;
        assert(sp.size() == ct_funs.size());

        int L = sp.size();

        std::vector<ce_vec> ce_levels(L);

        for (int l = 0; l < L; l++)
        {
            int I = sp[l].k_rho_vals.size();
            std::vector<cmplx> y(I);
            for (int i = 0; i < I; i++)
            {
                y[i] = gf.f(sp[l].k_rho_vals[i]) / (1.0i * sp[l].k_z_vals[i]);
                for (int j = 0; j < l - 1; j++)
                    y[i] -= eval_fun(ce_levels[j], sp[l].k_z_vals[i]);
            }
            auto ce_gpof = gpof::gpof(y, sp[l].d_t);
            ce_levels[l] = ct_funs[l](ce_gpof);
        }
        return ce_levels;
    }

    std::vector<CmplxImg> get_images(const std::vector<ce_vec> &ce_levels,
                                     real rho)
    {
        std::vector<CmplxImg> img_all;
        for (const auto &ce_level : ce_levels)
        {
            for (const auto &ce : ce_level)
            {
                auto r = calc_r(rho, ce.back());
                img_all.emplace_back(ce.front(), r);
            }
        }

        return img_all;
    }

    cmplx get_spectral_gf(const std::vector<ce_vec> &ce_levels,
                          real rho,
                          real k_0)
    {
        using std::complex_literals::operator""i;

        auto img_all = get_images(ce_levels, rho);

        cmplx val = 0;
        for (const auto &img : img_all)
            val += img.amplitude * std::exp(-1.0i * k_0 * img.r) / img.r;

        return val;
    }

} // namespace mthesis::dcim::utils
