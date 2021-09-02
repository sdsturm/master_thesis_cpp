#ifndef MTHESIS_GEOMETRY_HPP
#define MTHESIS_GEOMETRY_HPP

#include <mthesis/definitions.hpp>

#include <armadillo>

namespace mthesis
{
    using VectorR3 = arma::vec3;
    using VectorC3 = arma::cx_vec3;

    struct LayeredMediumCoords
    {
        const VectorR3 R;
        const real rho;
        const real phi;
        const real z;
        const real z_;

        LayeredMediumCoords(const VectorR3 &r, const VectorR3 &r_)
            : R(r - r_),
              rho(std::sqrt(std::pow(R[0], 2) + std::pow(R[1], 2))),
              phi(std::atan2(R[1], R[0])),
              z(r[2]),
              z_(r_[2])
        {
        }
    };

} // namespace mthesis

#endif
