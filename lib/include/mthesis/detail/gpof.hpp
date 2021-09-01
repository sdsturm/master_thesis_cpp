#ifndef MTHESIS_GPOF_HPP
#define MTHESIS_GPOF_HPP

#include <mthesis/detail/definitions.hpp>

#include <vector>
#include <array>

namespace mthesis::gpof
{
    struct Params
    {
        double tol;
        int M;
        int M_max;
        int L;

        Params();
        void set_tol(double tol);
        void set_M(int M);
        void set_M_max(int M_max);
        void set_L(int L);
    };

    std::vector<CmplxExp> gpof(const std::vector<cmplx> &y,
                               real d_t,
                               const Params params = Params());

    std::vector<cmplx> reconstruct_signal(std::vector<CmplxExp> &ce,
                                          double d_t,
                                          unsigned N);

} // namespace mthesis::gpof

#endif
