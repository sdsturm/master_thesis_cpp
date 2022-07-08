#include <mthesis/definitions.hpp>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <iostream>

using namespace mthesis;

int main()
{
    using std::complex_literals::operator""i;

    real t = 0;
    auto f = [&](cmplx z)
    {
        return exp(1.0i * t * z) / (pow(z, 2) + 1.0);
    };

    auto gamma = [](real t) { return 1.0i + exp(2.0i * M_PI * t); };

    auto gamma_ = [](real t) { return 2.0i * M_PI * exp(2.0i * M_PI * t); };

    auto integrand = [&](real t) { return f(gamma(t)) * gamma_(t); };

    using boost::math::quadrature::gauss_kronrod;

    cmplx num = gauss_kronrod<double, 15>::integrate(integrand, 0, 1);

    cmplx ref = exp(-t) / 2.0i;

    std::cout << "ref = " << ref << "\n";
    std::cout << "num = " << num / (2.0i * M_PI) << "\n";

    return 0;
}
