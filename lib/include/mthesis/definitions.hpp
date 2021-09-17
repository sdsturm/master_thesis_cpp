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

enum class RiemannSheet {I, II, III, IV};

using SpectralGF = std::function<cmplx(real z, real z_, cmplx k_rho)>;

struct CmplxExp
{
    CmplxExp(cmplx amplitude, cmplx exponent) : amp(amplitude), exp(exponent) {}

    const cmplx amp;
    const cmplx exp;
};

struct LMCoords
{
    LMCoords(const VectorR3 &r, const VectorR3 &r_);

    const VectorR3 R;
    const real rho;
    const real phi;
    const real z;
    const real z_;
};

cmplx cmplx_length(const VectorC3 &r);

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
