#ifndef MTHESIS_DCIM_UTILS_HPP
#define MTHESIS_DCIM_UTILS_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/spectral_gf.hpp>

#include <functional>

namespace mthesis::dcim::utils
{
    using ce_vec = std::vector<CmplxExp>;
    using ct_fun = std::function<ce_vec(const ce_vec &)>;

    struct CmplxImg
    {
        const cmplx amplitude;
        const cmplx r;

        CmplxImg(cmplx amplitude, cmplx r) : amplitude(amplitude), r(r) {}
    };

    std::vector<cmplx> get_k_rho_vals(const std::vector<cmplx> &k_z_vals,
                                      real k_0);

    struct SamplingPath
    {
        const std::vector<cmplx> k_z_vals;
        const std::vector<cmplx> k_rho_vals;
        const real d_t;
        const unsigned N;

        SamplingPath(real k_0, const std::vector<cmplx> &k_z_vals, real d_t);
    };

    real get_k_0(const LayeredMedium &lm);

    real find_k_max(const LayeredMedium &lm);

    real get_d_t(real T, int N);

    cmplx calc_r(real rho, cmplx alpha);

    cmplx eval_fun(const ce_vec &ce, cmplx k_z);

    std::vector<ce_vec> algo(const si::SpectralGF &gf,
                             const std::vector<SamplingPath> &sp,
                             const std::vector<ct_fun> &ct_funs);

    std::vector<CmplxImg> get_images(const std::vector<ce_vec> &ce_levels,
                                     real rho);

    cmplx get_spectral_gf(const std::vector<ce_vec> &ce_levels,
                          real rho,
                          real k_0);

} // namespace mthesis::dcim::utils

#endif