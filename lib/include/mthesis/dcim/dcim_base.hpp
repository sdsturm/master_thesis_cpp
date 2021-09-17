#ifndef MTHESIS_DCIM_BASE_HPP
#define MTHESIS_DCIM_BASE_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/si/sommerfeld_integral.hpp>

#include <functional>

namespace mthesis::dcim {

using CeVec = std::vector<CmplxExp>;
using CtFun = std::function<CeVec(const CeVec &)>;


struct SamplingPath
{
    SamplingPath(real k_0, const std::vector<cmplx> &k_z_vals, real d_t);

    const std::vector<cmplx> k_z_vals;
    const std::vector<cmplx> k_rho_vals;
    const real d_t;

private:
    static std::vector<cmplx> get_k_rho_vals(const std::vector<cmplx> &k_z_vals,
                                             real k_0);
};


struct CmplxImg
{
    CmplxImg(cmplx amplitude, cmplx r);

    const cmplx amp;
    const cmplx r;
};


class DCIM
{
public:
    DCIM(const SommerfeldIntegral &si);

    // Make this a an abstract base class.
    virtual ~DCIM() = 0;

    std::vector<CeVec> get_exponentials(real z, real z_) const;

    std::vector<CmplxImg> get_images(const std::vector<CeVec> ce_vecs_levels,
                                     real rho) const;

    std::vector<CmplxImg> get_images(real z, real z_, real rho) const;

    cmplx get_spatial_gf(const std::vector<CeVec> ce_vecs_levels,
                         real rho) const;

    cmplx get_spatial_gf(real z, real z_, real rho) const;

    cmplx get_spatial_gf(const VectorR3 &r, const VectorR3 &r_) const;

protected:
    const SommerfeldIntegral &si;
    real k_0;
    real k_max;
    std::vector<SamplingPath> sampling_paths;
    std::vector<CtFun> ct_funs;
};


// Hide implementation details in nested namespace utils.
namespace utils {

real get_d_t(real T, int N);

cmplx calc_r(real rho, cmplx alpha);

cmplx eval_fun(const CeVec &ce, cmplx k_z);


std::vector<cmplx> get_y(const std::function<cmplx(cmplx)> &G,
                         const std::vector<CeVec> &ce_levels,
                         const std::vector<SamplingPath> &sampling_paths,
                         int level);

} // namespace utils

} // namespace mthesis::dcim

#endif // MTHESIS_DCIM_BASE_HPP
