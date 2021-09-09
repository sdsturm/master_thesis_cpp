#ifndef MTHESIS_DCIM_TWO_LEVEL_HPP
#define MTHESIS_DCIM_TWO_LEVEL_HPP

#include <mthesis/dcim/dcim_base.hpp>

namespace mthesis::dcim {

class TwoLevel : public DCIM
{
public:
    TwoLevel(const SommerfeldIntegral &si);

private:
SamplingPath calc_sp_2(real T_2, int N_2) const;

SamplingPath calc_sp_1(real T_1, real T_2, int N_1) const;

CeVec calc_coeffs_1(const CeVec &ce_in, real T_2) const;

CeVec calc_coeffs_2(const CeVec &ce_in, real T_2) const;
};

} // namespace mthesis::dcim

#endif // MTHESIS_DCIM_TWO_LEVEL_HPP
