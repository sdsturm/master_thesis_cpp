#ifndef MTHESIS_SI_PE_PARAMS_HPP
#define MTHESIS_SI_PE_PARAMS_HPP

#include <cassert>

namespace mthesis::si::pe {

struct Params
{
    unsigned max_intervals;
    double tol;

    Params() : max_intervals(15), tol(1e-8) {}

    void set_max_intervals(int max_intervals)
    {
        assert(max_intervals > 0);
        this->max_intervals = max_intervals;
    }

    void set_tol(double tol)
    {
        assert(tol > 0.0);
        this->tol = tol;
    }
};

} // namespace mthesis::si::pe

#endif
