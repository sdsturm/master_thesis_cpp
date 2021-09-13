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

    auto g_1 = [&](const VectorR3 &r, const VectorC3 &r_)
    {
        VectorC3 R = r - r_;
        cmplx R_cl = std::sqrt(arma::dot(R, R));
        if (R_cl.real() < 0.0) {
            R_cl *= -1;
        }

        return exp(-1.0i * fd.k_0 * R_cl) / R_cl;
    };

    auto g_2 = [&](const VectorR3 &r, const VectorC3 &r_)
    {
        VectorC3 R = r - r_;
        cmplx R_cl = std::sqrt(arma::dot(R, R));
        if (R_cl.imag() > 0.0) {
            R_cl *= -1;
        }

        return exp(-1.0i * fd.k_0 * R_cl) / R_cl;
    };

    unsigned N_pts = 60;
    arma::vec x_vals = arma::linspace(-5, 5, N_pts) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-5, 5, N_pts) * fd.lambda_0;

    VectorC3 r_1_ = {0, 0, -0.01i * fd.lambda_0};
    VectorC3 r_2_ = {0, 0, -0.1i * fd.lambda_0};
    VectorC3 r_3_ = {0, 0, -1.0i * fd.lambda_0};

    printf("x_by_lambda_0 z_by_lambda_0 g_1_1 g_1_2 g_1_3 g_2_1 g_2_2 g_2_3\n");
    for (const auto &x : x_vals) {
        for (const auto &z : z_vals) {
            VectorR3 r = {x, 0, z};
            cmplx g_1_1 = g_1(r, r_1_);
            cmplx g_1_2 = g_1(r, r_2_);
            cmplx g_1_3 = g_1(r, r_3_);
            cmplx g_2_1 = g_2(r, r_1_);
            cmplx g_2_2 = g_2(r, r_2_);
            cmplx g_2_3 = g_2(r, r_3_);

            printf("%.8f %.8f %.8e %.8e %.8e %.8e %.8e %.8e\n",
                   x / fd.lambda_0,
                   z / fd.lambda_0,
                   g_1_1.real(),
                   g_1_2.real(),
                   g_1_3.real(),
                   g_2_1.real(),
                   g_2_2.real(),
                   g_2_3.real()
                   );
        }
        printf("\n");
    }

    return 0;
}
