#ifndef MTHESIS_FMM_GEGENBAUER_HPP
#define MTHESIS_FMM_GEGENBAUER_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>
#include <mthesis/fmm/helpers.hpp>

namespace mthesis::fmm {

cmplx compute_gegenbauer_pw(const Params &params,
                            const VectorR3 &r,
                            const VectorR3 &r_,
                            unsigned L);

// Overload for complex source point.
cmplx compute_gegenbauer_pw(const Params &params,
                            const VectorR3 &r,
                            const VectorC3 &r_,
                            unsigned L);

} // namespace mthesis::fmm

#endif // MTHESIS_FMM_GEGENBAUER_HPP
