#ifndef MTHESIS_DEFINIITONS_HPP
#define MTHESIS_DEFINIITONS_HPP

#include <armadillo>

#include <complex>
#include <array>

namespace mthesis
{
    using real = double;
    using cmplx = std::complex<real>;
    using VectorR3 = arma::vec3;
    using VectorC3 = arma::cx_vec3;

    using CmplxExp = std::array<cmplx, 2>;

    enum class EmMode
    {
        TM,
        TE
    };

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
}

#endif // namespace mthesis
