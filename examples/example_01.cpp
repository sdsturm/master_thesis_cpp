#include <mthesis/mthesis.hpp>

#include <iostream>

using namespace mthesis;

int main()
{
    auto fd = FrequencyDomain(1e9);
    auto eps_r = cmplx_permittivity(fd, 3, 1e-3);
    auto ground = Medium(fd, eps_r, 1.0);
    auto hs = HalfSpace(fd, ground);

    std::cout << fd.lambda_0 << "\n";
    std::cout << ground.k << "\n";

    return 0;
}
