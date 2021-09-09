#include <mthesis.hpp>

#include <armadillo>
#include <gsl/gsl_const_mksa.h>

#include <cstdio>

using namespace mthesis;

int main()
{
    using std::complex_literals::operator""i;

    FrequencyDomain fd(GSL_CONST_MKSA_SPEED_OF_LIGHT);
    Vacuum medium(fd);

    unsigned N_pts = 80;
    arma::vec x_vals = arma::linspace(-5, 5, N_pts) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-5, 5, N_pts) * fd.lambda_0;

    VectorC3 r_ = {(0.0 + 0.1i) * fd.lambda_0, 0, 0};

    printf("x_by_lambda_0 z_by_lambda_0 g_re g_im\n");
    for (const auto &x : x_vals) {
        for (const auto &z : z_vals) {
            VectorR3 r = {x, 0, z};
            cmplx val = gf::scalar::free_space::G_0(medium, r, r_);
            printf("%.8f %.8f %.8e %.8e\n",
                   x / fd.lambda_0,
                   z / fd.lambda_0,
                   val.real(),
                   val.imag());
        }
        printf("\n");
    }

    return 0;
}
