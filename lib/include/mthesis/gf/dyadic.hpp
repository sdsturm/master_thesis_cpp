#ifndef MTHESIS_GF_DYADIC_HPP
#define MTHESIS_GF_DYADIC_HPP

#include <mthesis/definitions.hpp>
#include <mthesis/solution_domain.hpp>

namespace mthesis::gf::dyadic
{

namespace free_space
{

DyadC3 G_EJ(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_EM(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HJ(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HM(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

// Hide implementation details in nested namespace utils.
namespace utils {

DyadC3 G_e0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_m0(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

} // namespace utils

} // namespace free_space

namespace layered_media
{

DyadC3 G_EJ(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_EM(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HJ(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

DyadC3 G_HM(const Medium &medium, const VectorR3 &r, const VectorR3 &r_);

// Hide implementation details in nested namespace utils.
namespace utils {

// TODO

} // namespace utils

} // namespace layered_media

} // namespace mthesis::gf::dyadic

#endif // MTHESIS_GF_DYADIC_HPP
