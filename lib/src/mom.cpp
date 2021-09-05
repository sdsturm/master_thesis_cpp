#include <mthesis/mom.hpp>

#include <cassert>

namespace mthesis::mom {

MoM::MoM(const FrequencyDomain &fd,
         const std::vector<VectorR3> &src_pts,
         const std::vector<VectorR3> &obs_pts)
    : A(obs_pts.size(), src_pts.size())
{
    using std::complex_literals::operator""i;

    for (size_t m = 0; m < obs_pts.size(); m++) {
        for (size_t m_ = 0; m_ < src_pts.size(); m_++) {
            auto R = arma::norm(obs_pts[m] - src_pts[m_]);
            A(m, m_) = std::exp(1.0i * fd.k_0 * R) / R;
        }
    }
}

std::vector<cmplx> MoM::calc_product(const std::vector<cmplx> &I) const
{
    assert(I.size() == A.n_cols);

    const arma::cx_vec I_map(I);

    std::vector<cmplx> V(A.n_rows);
    arma::cx_vec V_map(V.data(), V.size(), false, true);

    V_map = A * I_map;

    return V;
}

} // namespace mthesis::mom
