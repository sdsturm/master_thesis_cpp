#include "./misc/print_dgf_half_space.hpp"

using namespace mthesis;

int main()
{
    VectorR3 J = {0, 0, 1}; // VED current.
    print_dgf_half_space(J);

    return 0;
}
