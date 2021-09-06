#include <mthesis.hpp>

#include <boost/math/special_functions/legendre.hpp>
#include <armadillo>

#include <cstdio>

int main()
{
#if 1
    arma::vec x_vals = arma::linspace(-1, 1, 51);

    unsigned nu_max = 500;
    for (unsigned nu = 0; nu <= nu_max; nu++) {
        double err_max = 0;
        for (const auto &x : x_vals) {
            auto ref = boost::math::legendre_p(nu, x);
            auto num = mthesis::legendre_p_recurrence(nu, x);
            auto err = std::abs(num - ref);
            if (err > err_max) {
                err_max = err;
            }
        }
        printf("L = %3d    max_err = %.4e\n", nu, err_max);
    }
#else
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
#endif

    return 0;
}
