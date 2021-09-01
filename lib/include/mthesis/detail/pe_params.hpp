#ifndef MTHESIS_PE_PARAMS_HPP
#define MTHESIS_PE_PARAMS_HPP

namespace mthesis
{
    struct PEParams
    {
        const unsigned max_intervals;
        const double tol;
        PEParams(unsigned max_intervals = 15, double tol = 1e-8)
            : max_intervals(max_intervals), tol(tol)
        {
        }
    };
    
} // namespace mthesis


#endif
