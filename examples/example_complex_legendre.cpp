#include <mthesis.hpp>

#include <boost/math/special_functions/legendre.hpp>
#include <armadillo>

#include <cstdio>

int main()
{
    arma::vec x_vals = arma::linspace(-1, 1, 61);
    arma::vec y_vals = arma::linspace(-1, 1, 61);

    size_t N = x_vals.n_elem * y_vals.n_elem;

    unsigned nu_max = 5;

    // Print header.
    printf("z_re z_im ");
    for (unsigned nu = 0; nu <= nu_max; nu++) {
        printf("nu_%d_abs nu_%d_arg ", nu, nu);
    }
    printf("\n");

    size_t n = 0;
    for (const auto &x : x_vals) {
        for (const auto &y : y_vals) {
            mthesis::cmplx z(x, y);
            printf("%.6f %.6f ", x, y);

            for (unsigned nu = 0; nu <= nu_max; nu++) {
                auto ans = mthesis::legendre_p_recurrence(nu, z);
                printf("%.6e %.6e ", std::abs(ans), std::arg(ans));
            }
            printf("\n");
        }
        printf("\n");
    }

    return 0;
}
