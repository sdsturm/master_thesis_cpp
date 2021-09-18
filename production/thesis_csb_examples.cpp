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

    auto g1 = [&](const VectorR3 &r, const VectorC3 &r_)
    {
        VectorC3 rdiff = r - r_;

        cmplx R;
        for (const auto &comp: rdiff) {
            R += pow(comp, 2);
        }
        R = sqrt(R);

        // Select branch with positive real part.
        if (R.real() < 0) {
            R *= -1;
        }

        return exp(-1.0i * fd.k_0 * R) / (4.0 * M_PI * R);
    };

    auto g2 = [&](const VectorR3 &r, const VectorC3 &r_)
    {
        VectorC3 rdiff = r - r_;

        cmplx R;
        for (const auto &comp: rdiff) {
            R += pow(comp, 2);
        }
        R = sqrt(R);

        // Select branch with negative imaginary part.
        if (R.imag() > 0) {
            R *= -1;
        }

        return exp(-1.0i * fd.k_0 * R) / (4.0 * M_PI * R);
    };

    unsigned N_pts = 60;
    arma::vec x_vals = arma::linspace(-5, 5, N_pts) * fd.lambda_0;
    arma::vec z_vals = arma::linspace(-5, 5, N_pts) * fd.lambda_0;

    VectorC3 r_1_ = {0, 0, -0.1i * fd.lambda_0};
    VectorC3 r_2_ = {0, 0, -0.5i * fd.lambda_0};
    VectorC3 r_3_ = {0, 0, -1.0i * fd.lambda_0};

    printf("x_by_lambda_0 z_by_lambda_0 g1_1 g1_2 g1_3 g2_1 g2_2 g2_3\n");
    for (const auto &x : x_vals) {
        for (const auto &z : z_vals) {
            VectorR3 r = {x, 0, z};
            cmplx g1_1 = g1(r, r_1_);
            cmplx g1_2 = g1(r, r_2_);
            cmplx g1_3 = g1(r, r_3_);
            cmplx g2_1 = g2(r, r_1_);
            cmplx g2_2 = g2(r, r_2_);
            cmplx g2_3 = g2(r, r_3_);

            printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                   x / fd.lambda_0,
                   z / fd.lambda_0,
                   g1_1.real(),
                   g1_2.real(),
                   g1_3.real(),
                   g2_1.real(),
                   g2_2.real(),
                   g2_3.real()
                   );
        }
        printf("\n");
    }

    return 0;
}
