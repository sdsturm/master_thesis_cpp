#ifndef MTHESIS_DCIM_DCIM_HPP
#define MTHESIS_DCIM_DCIM_HPP

#include <mthesis/detail/dcim_utils.hpp>

namespace mthesis::dcim::threelevelv2
{
    real calc_T_03(real k_0, real k_rho_max_3);

    real calc_T_02(real k_0, real k_rho_max_2, real T_03);

    real calc_T_01(real k_0, real k_rho_max_1, real T_02, real T_03);

    utils::SamplingPath calc_C_3(real k_0, real T_03, int N);

    utils::SamplingPath calc_C_2(real k_0, real T_02, real T_03, int N);

    utils::SamplingPath calc_C_1(real k_0, real T_01, real T_02, real T_03, int N);

    utils::ce_vec calc_coeffs_1(const utils::ce_vec &ce_in,
                                real k_0,
                                real T_02,
                                real T_03);

    utils::ce_vec calc_coeffs_2(const utils::ce_vec &ce_in,
                                real k_0,
                                real T_02,
                                real T_03);

    utils::ce_vec calc_coeffs_3(const utils::ce_vec &ce_in,
                                real k_0,
                                real T_03);

    std::vector<utils::ce_vec> three_level_v2(const si::SpectralGF &gf);

} // namespace mthesis::dcim::threelevelv2

#endif
