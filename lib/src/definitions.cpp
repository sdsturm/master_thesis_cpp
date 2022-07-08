#include <mthesis/definitions.hpp>

namespace mthesis {

LMCoords::LMCoords(const VectorR3 &r, const VectorR3 &r_)
    : R(r - r_),
      rho(std::sqrt(std::pow(R[0], 2) + std::pow(R[1], 2))),
      phi(std::atan2(R[1], R[0])),
      z(r[2]),
      z_(r_[2])
{}

cmplx cmplx_length(const VectorC3 &r)
{
    cmplx val;
    for (const auto &component : r) {
        val += std::pow(component, 2.0);
    }
    val = std::sqrt(val);

    // Choose branch with non-negative real part, see Hansen2013.
    if (val.real() < 0.0) {
        val *= -1.0;
    }

    return val;
}

} // namespace mthesis
