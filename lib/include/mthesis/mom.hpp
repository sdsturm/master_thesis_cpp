#ifndef MTHESIS_MOM_HPP
#define MTHESIS_MOM_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

namespace mthesis::mom {

class MoM
{
public:
    MoM(const FrequencyDomain &fd,
        const std::vector<VectorR3> &src_pts,
        const std::vector<VectorR3> &obs_pts);

    std::vector<cmplx> calc_product(const std::vector<cmplx> &I) const;

private:
    arma::cx_mat A;
};

} // namespace mthesis::mom

#endif // MTHESIS_MOM_HPP
