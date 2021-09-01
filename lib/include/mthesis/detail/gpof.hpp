#ifndef MTHESIS_GPOF_HPP
#define MTHESIS_GPOF_HPP

#include <mthesis/detail/definitions.hpp>

#include <vector>

namespace mthesis
{
    struct ComplexExponential
    {
        cmplx amplitude;
        cmplx exponent;

        ComplexExponential(cmplx amplitude, cmplx exponent)
            : amplitude(amplitude), exponent(exponent)
        {
        }
    };

    struct GPOFParams
    {
        double tol;
        int M;
        int M_max;
        int L;

        GPOFParams(double tol = 1e-6, int M = -1, int M_max = 10, int L = -1)
            : tol(tol), M(M), M_max(M_max), L(L)
        {
        }
    };

    std::vector<ComplexExponential> gpof(const std::vector<cmplx> &y,
                                         real d_t,
                                         const GPOFParams params = GPOFParams());

    std::vector<cmplx> reconstruct_signal(std::vector<ComplexExponential> &ce,
                                          double d_t,
                                          unsigned N);

} // namespace mthesis

#endif
