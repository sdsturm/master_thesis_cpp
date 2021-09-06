#ifndef MTHESIS_GF_DYADIC_HPP
#define MTHESIS_GF_DYADIC_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

namespace mthesis::gf::dyadic {

namespace free_space {

DyadC3 G_EJ(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_EM(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HJ(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HM(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

// Hide implementation details in nested namespace utils.
namespace utils {

DyadC3 G_e0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_m0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

} // namespace utils

} // namespace free_space

namespace layered_media
{

DyadC3 G_EJ(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HM(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HJ(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_EM(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_);

// Hide implementation details in nested namespace utils.
namespace utils {

real get_Z_0();

DyadC3 G_EJ_HM_common(const std::vector<cmplx> &si_vals,
                      const LayeredMediumCoords &coords,
                      cmplx f,
                      cmplx f_);

cmplx calc_electric_factor(const LayeredMedium &lm, real z);

cmplx calc_magnetic_factor(const LayeredMedium &lm, real z);

std::vector<cmplx> calc_G_EJ_si_vals(const LayeredMedium &lm,
                                     const LayeredMediumCoords &coords);

std::vector<cmplx> calc_G_HM_si_vals(const LayeredMedium &lm,
                                     const LayeredMediumCoords &coords);

std::vector<cmplx> calc_G_HJ_si_vals(const LayeredMedium &lm,
                                     const LayeredMediumCoords &coords);

cmplx G_EJ_HM_xx(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords);

cmplx G_EJ_HM_yx(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords);

cmplx G_EJ_HM_xz(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords,
                 cmplx f_);

cmplx G_EJ_HM_yy(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords);

cmplx G_EJ_HM_yz(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords,
                 cmplx f_);

cmplx G_EJ_HM_zx(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords,
                 cmplx f);

cmplx G_EJ_HM_zy(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords,
                 cmplx f);

cmplx G_EJ_HM_zz(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords,
                 cmplx f,
                 cmplx f_);

cmplx G_HJ_xx(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords);

cmplx G_HJ_xy(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords);

cmplx G_HJ_xz(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords,
              cmplx factor_e_);

cmplx G_HJ_yx(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords);

cmplx G_HJ_yz(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords,
              cmplx factor_e_);

cmplx G_HJ_zx(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords,
              cmplx factor_h);

cmplx G_HJ_zy(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords,
              cmplx factor_h);

} // namespace utils

} // namespace layered_media

} // namespace mthesis::gf::dyadic

#endif // MTHESIS_GF_DYADIC_HPP
