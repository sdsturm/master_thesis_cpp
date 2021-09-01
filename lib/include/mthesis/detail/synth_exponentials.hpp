#ifndef MTHESIS_SYNTH_EXPONENTIALS_HPP
#define MTHESIS_SYNTH_EXPONENTIALS_HPP

#include <mthesis/detail/gpof.hpp>

#include <random>

namespace mthesis
{
    class SynthExponentials
    {
    private:
        unsigned N;
        double d_t;
        std::vector<double> t;
        std::vector<cmplx> y;
        std::vector<ComplexExponential> ce;

    public:
        SynthExponentials(unsigned M, unsigned N, double oversampling)
            : N(N), t(N), y(N)
        {
            constexpr double f_max = 30;
            double f_s = f_max * oversampling;
            double d_t = 1.0 / f_s;

            for (unsigned n = 0; n < N; n++)
                this->y[n] = n * d_t;

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
                theta[m] =dist_theta(gen);
            }
            f.front() = 0.0;

            using std::complex_literals::operator""i;
            for (size_t m = 0; m < M; m++)
            {
                cmplx amplitude = A[m] * std::exp(1.0i * theta[m]);
                cmplx exponent = alpha[m] + 2.0i * M_PI * f[m];
                this->ce.emplace_back(amplitude, exponent);
            }

            this->y = reconstruct_signal(this->ce, this->d_t, this->N);
        }

        unsigned get_N() const { return N; }
        unsigned get_M() const { return ce.size(); }
        double get_d_t() const { return d_t; }
        std::vector<double> get_t() const { return t; }
        std::vector<cmplx> get_y() const { return y; }
    };

} // namespace mthesis

#endif
