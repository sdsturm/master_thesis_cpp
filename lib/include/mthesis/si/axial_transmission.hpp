#ifndef MTHESIS_SI_AXIAL_TRANSMISSION_HPP
#define MTHESIS_SI_AXIAL_TRANSMISSION_HPP

#include <mthesis/si/sommerfeld_integral.hpp>

namespace mthesis::si::axial_transmission {

cmplx eval_si_along_sip(const SommerfeldIntegral &si,
                        real rho,
                        real z,
                        real z_,
                        pe::Params pe_params = pe::Params());

cmplx eval_si_along_sip(const SommerfeldIntegral &si,
                        const VectorR3 &r,
                        const VectorR3 &r_,
                        pe::Params pe_params = pe::Params());

namespace utils {

real calc_pe_start(const SommerfeldIntegral &si);

real calc_indention(real rho, real a);

cmplx eval_head_ellipsis(const SommerfeldIntegral &si,
                         real rho,
                         real z,
                         real z_,
                         real a);

cmplx eval_tail_on_sip(const SommerfeldIntegral &si,
                       real rho,
                       real z,
                       real z_,
                       real a,
                       const pe::Params &pe_params);

} // namespace utils

} // namespace mthesis::si::axial_transmission

#endif // MTHESIS_SI_AXIAL_TRANSMISSION_HPP
