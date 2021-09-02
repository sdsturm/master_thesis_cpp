#ifndef MTHESIS_SYNTH_EXPONENTIALS_HPP
#define MTHESIS_SYNTH_EXPONENTIALS_HPP

#include <mthesis/gpof.hpp>

#include <random>

namespace mthesis::gpof
{
    struct SynthExponentials
    {
        unsigned N;
        double d_t;
        std::vector<double> t;
        std::vector<cmplx> y;
        std::vector<CmplxExp> ce;

        SynthExponentials(unsigned M, unsigned N, double oversampling)
            : N(N), t(N), y(N)
        {
            constexpr double f_max = 30;
            double f_s = f_max * oversampling;
            d_t = 1.0 / f_s;

            for (unsigned n = 0; n < N; n++)
            {
                t[n] = n * d_t;
            }

            std::random_device dev;
            std::mt19937 gen(dev());
            std::uniform_real_distribution<double> dist_A(1.0, 9.0);
            std::uniform_real_distribution<double> dist_alpha(-4.0, 0.0);
            std::uniform_real_distribution<double> dist_f(0.0, f_max);
            std::uniform_real_distribution<double> dist_theta(-M_PI, M_PI);

            std::vector<double> A(M), alpha(M), f(M), theta(M);
            for (unsigned m = 0; m < M; m++)
            {
                A[m] = dist_A(gen);
                alpha[m] = dist_alpha(gen);
                f[m] = dist_f(gen);
                theta[m] = dist_theta(gen);
            }
            f.front() = 0.0;

            using std::complex_literals::operator""i;
            for (size_t m = 0; m < M; m++)
            {
                cmplx b = A[m] * std::exp(1.0i * theta[m]);
                cmplx beta = alpha[m] + 2.0i * M_PI * f[m];
                ce.push_back(CmplxExp{b, beta});
            }

            y = reconstruct_signal(ce, d_t, N);
        }
    };

} // namespace mthesis::gpof

#endif
