#ifndef MTHESIS_FMM_FREE_SPACE_HPP
#define MTHESIS_FMM_FREE_SPACE_HPP

#include <mthesis/fmm/helpers.hpp>

namespace mthesis::fmm::freespace {

unsigned calc_L(const Params &params);

f_of_k_hat calc_ff(const FrequencyDomain &fd,
                   const EwaldSphere &es,
                   const Group &g,
                   const VectorR3 &r);

std::vector<f_of_k_hat> calc_all_ff(const FrequencyDomain &fd,
                                    const EwaldSphere &es,
                                    const std::vector<Group> &groups,
                                    const std::vector<VectorR3> &pts);

f_of_k_hat calc_top(unsigned L,
                    const FrequencyDomain &fd,
                    const EwaldSphere &es,
                    const Group &src_group,
                    const Group &obs_group);

std::vector<f_of_k_hat> calc_all_top(unsigned L,
                                     const FrequencyDomain &fd,
                                     const EwaldSphere &es,
                                     const std::vector<Group> &src_groups,
                                     const std::vector<Group> &obs_groups);

struct FMM
{
    const FrequencyDomain &fd;
    const real w;
    const std::vector<VectorR3> &src_pts;
    const std::vector<VectorR3> &obs_pts;
    const std::vector<Group> src_groups;
    const std::vector<Group> obs_groups;
    const unsigned L;
    const EwaldSphere es;
    const std::vector<f_of_k_hat> src_ff_all;
    const std::vector<f_of_k_hat> obs_ff_all;
    const std::vector<f_of_k_hat> top_all;

    FMM(const Params &params,
                 const std::vector<VectorR3> &src_pts,
                 const std::vector<VectorR3> &obs_pts);

    std::vector<cmplx> calc_product(const std::vector<cmplx> &I) const;
};

} // namespace mthesis::fmm::freespace

#endif // MTHESIS_FMM_FREE_SPACE_HPP
