#ifndef MTHESIS_DCIM_THREE_LEVEL_V3_HPP
#define MTHESIS_DCIM_THREE_LEVEL_V3_HPP

#include <mthesis/dcim/dcim_base.hpp>

namespace mthesis::dcim {

class ThreeLevelV3 : public DCIM
{
public:
    ThreeLevelV3(const SommerfeldIntegral &si);

private:
real calc_T_2(real k_rho_max_2) const;

real calc_T_1(real k_rho_max_1, real T_2) const;

SamplingPath calc_sp_3(real T_3, int N_3) const;

SamplingPath calc_sp_2(real T_2, int N_2) const;

SamplingPath calc_sp_1(real T_1, real T_2, int N_1) const;

CeVec calc_coeffs_1(const CeVec &ce_in, real T_2) const;

CeVec calc_coeffs_2(const CeVec &ce_in) const;

CeVec calc_coeffs_3(const CeVec &ce_in, real T_3) const;
};

} // namespace mthesis::dcim


#endif // MTHESIS_DCIM_THREE_LEVEL_V3_HPP
