#include "./misc/print_dgf_half_space.hpp"

using namespace mthesis;

int main()
{
    VectorR3 J = {1, 0, 0}; // VED current.
    print_dgf_half_space(J);

    return 0;
}
