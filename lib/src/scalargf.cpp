#include <mthesis/scalargf.hpp>

namespace mthesis::scalargf {

namespace freespace {

cmplx G_0(const Medium &m, const VectorR3 &r, const VectorR3 &r_)

} // namespace freespace

namespace layeredmedia {

cmplx generic_spectral(const LayeredMedium &lm,
                          real z,
                          real z_,
                          cmplx k_rho,
                          EmMode mode,
                          bool direct_term)

SommerfeldIntegral get_sommerfeld_integral(const LayeredMedium &lm,
                                           real nu,
                                           EmMode mode,
                                           bool direct_term,
                                           SiParams si_params)

} // namespace layeredmedia

} // namespace mthesis::scalargf
