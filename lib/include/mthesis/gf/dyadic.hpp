#ifndef MTHESIS_GF_DYADIC_HPP
#define MTHESIS_GF_DYADIC_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

namespace mthesis::gf::dyadic {

// *****************************************************************************
//                   Free-space dyadic Green's functions
// *****************************************************************************
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


// *****************************************************************************
//                  Layered media dyadic Green's function
// *****************************************************************************

namespace layered_media
{

DyadC3 G_EJ(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HM(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HJ(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_EM(const LayeredMedium &lm, const VectorR3 &r, const VectorR3 &r_);

// Hide implementation details in nested namespace utils.
namespace utils {

real get_eta_0();

cmplx calc_electric_factor(const LayeredMedium &lm, real z);

cmplx calc_magnetic_factor(const LayeredMedium &lm, real z);

// Note: G_EJ and G_HM follow the same schme and are duals of each other.
std::vector<cmplx> calc_G_EJ_si_vals(const LayeredMedium &lm,
                                     const LayeredMediumCoords &coords);

std::vector<cmplx> calc_G_HM_si_vals(const LayeredMedium &lm,
                                     const LayeredMediumCoords &coords);

cmplx G_EJ_HM_xx(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords);

// Note: E_EJ_HM_xy is equal to G_EJ_HM_yx.

cmplx G_EJ_HM_xz(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords,
                 cmplx f_);

cmplx G_EJ_HM_yx(const std::vector<cmplx> &si_vals,
                 const LayeredMediumCoords &coords);

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

DyadC3 G_EJ_HM_common(const std::vector<cmplx> &si_vals,
                      const LayeredMediumCoords &coords,
                      cmplx f,
                      cmplx f_);

// Note: G_HJ and G_EM can be obtained from each other by reciprocity.
std::vector<cmplx> calc_G_HJ_si_vals(const LayeredMedium &lm,
                                     const LayeredMediumCoords &coords);

cmplx G_HJ_xx(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords);

cmplx G_HJ_xy(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords);

cmplx G_HJ_xz(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords,
              cmplx factor_e_);

cmplx G_HJ_yx(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords);

// Note: G_HJ_yy is -G_HJ_xx.

cmplx G_HJ_yz(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords,
              cmplx factor_e_);

cmplx G_HJ_zx(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords,
              cmplx factor_h);

cmplx G_HJ_zy(const std::vector<cmplx> &si_vals,
              const LayeredMediumCoords &coords,
              cmplx factor_h);

// Note: G_HJ_zz is always zero.

} // namespace utils

} // namespace layered_media

} // namespace mthesis::gf::dyadic

#endif // MTHESIS_GF_DYADIC_HPP
