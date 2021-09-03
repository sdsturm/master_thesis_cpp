#ifndef MTHESIS_FMM_HPP
#define MTHESIS_FMM_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

#include <armadillo>

namespace mthesis::fmm {

// Definition: FMM group indexing scheme.
//              |       |       |       |       |       |       |
// group index: |  -2   |  -1   |   0   |   1   |   2   |   3   |
//              |       |       |       |       |       |       |
// -----|-------|-------|-------|-------|-------|-------|-------|-----> x, y, z
//     -3w     -2w     -1w      0      +1w      2w      3w     4w

// Convention: a point (x, y, z) lies in the i-th group of the corresponding
// spatial direction if i * w <= x, y, z < (i+1) * w.

struct Params
{
    const FrequencyDomain &fd;
    const real w;

    Params(const FrequencyDomain &fd, real w);
};

using multiindex = arma::ivec3;

VectorR3 group_center(const Params &params, const multiindex &mi);

multiindex identify_group(const Params &params, const VectorR3 &r);

std::vector<VectorR3> rand_pts_in_group(const Params &params,
                                        const multiindex &mi,
                                        unsigned N);

void append_pts(std::vector<VectorR3> &pts,
                const std::vector<VectorR3> &new_pts);

struct Group
{
    const multiindex mi;
    const VectorR3 r_center;
    std::vector<size_t> pts_idx;

    Group(const Params &params, const multiindex &mi);

    void add_point_index(size_t index);
};

std::vector<Group> build_groups(const Params &params,
                                const std::vector<VectorR3> &pts);

class GroupSeparationError : public std::exception
{
private:
    const Group sg;
    const Group og;

public:
    GroupSeparationError(const Group &sg, const Group &og);
    void print_message() const;
};

std::vector<Group> build_groups(const Params &params,
                                const std::vector<VectorR3> &pts);

struct EwaldSpere
{
    std::vector<real> cos_theta_weights;
    real prefactor;
    std::vector<VectorR3> k_hat;

    EwaldSpere(unsigned L);
    cmplx integrate(const std::vector<cmplx> &f_of_k_hat) const;
};

unsigned calc_L(const Params &params);

std::vector<cmplx> calc_ff(const FrequencyDomain &fd,
                           const EwaldSpere &es,
                           const Group &g,
                           const VectorR3 &r);

std::vector<cmplx> calc_top(unsigned L,
                            const FrequencyDomain &fd,
                            const EwaldSpere &es,
                            const Group &src_group,
                            const Group &obs_group);

struct FreeSpaceFMM
{
    const FrequencyDomain &fd;
    const real w;
    const std::vector<VectorR3> &src_pts;
    const std::vector<VectorR3> &obs_pts;
    const std::vector<Group> src_groups;
    const std::vector<Group> obs_groups;
    const unsigned L;
    const EwaldSpere es;

    FreeSpaceFMM(const Params &params,
                 const std::vector<VectorR3> &src_pts,
                 const std::vector<VectorR3> &obs_pts);


};

} // namespace mthesis::fmm

#endif // MTHESIS_FMM_HPP
