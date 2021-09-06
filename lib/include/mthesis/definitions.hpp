#ifndef MTHESIS_DEFINIITONS_HPP
#define MTHESIS_DEFINIITONS_HPP

#include <armadillo>

#include <complex>
#include <cmath>

namespace mthesis {

using real = double;
using cmplx = std::complex<real>;

enum class EmMode{TM, TE};

using VectorR3 = arma::vec3;
using VectorC3 = arma::cx_vec3;
using DyadC3 = arma::cx_mat33;

struct CmplxExp
{
    const cmplx amp;
    const cmplx exp;
    CmplxExp(cmplx amplitude, cmplx exponent) : amp(amplitude), exp(exponent) {}
};

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

inline cmplx cmplx_length(const VectorC3 &r)
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

template <typename T1, typename T2>
real calc_rel_err(T1 num, T2 ref)
{
    return std::abs((num - ref) / ref);
}

template <typename T1, typename T2>
real calc_rel_err_db(T1 num, T2 ref)
{
    return 20.0 * std::log10(calc_rel_err(num, ref));
}

} // namespace mthesis

#endif
