#ifndef MTHESIS_WRITE_GENERIC_SGF_NEAR_REGION_HPP
#define MTHESIS_WRITE_GENERIC_SGF_NEAR_REGION_HPP

#include <mthesis.hpp>

#include <armadillo>

#include <filesystem>
#include <cstdio>
#include <cassert>

namespace mthesis {

void write_generic_sgf_near_region(const FrequencyDomain &fd,
                                   cmplx eps_r,
                                   const arma::vec &rho_vals_by_lambda_0,
                                   const arma::vec &z_vals_by_lambda_0,
                                   const std::filesystem::path &target_dir,
                                   const std::string &case_name)
{
    assert(arma::min(rho_vals_by_lambda_0) >= 0.0);

}

} // namespace mthesis

#endif // MTHESIS_WRITE_GENERIC_SGF_NEAR_REGION_HPP
