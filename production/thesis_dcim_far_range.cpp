#include "./misc/print_dcim_error.hpp"

#include <armadillo>

int main()
{
    arma::vec rho_vals_by_lambda_0 = arma::logspace(-3, 3, 60);
    arma::vec z_vals_by_lambda_0 = arma::logspace(-3, 3, 60);

    mthesis::print_dcim_error(rho_vals_by_lambda_0, z_vals_by_lambda_0);

    return 0;
}
