#include <mthesis/detail/definitions.hpp>
#include <mthesis/detail/gpof.hpp>
#include <mthesis/detail/synth_exponentials.hpp>

#include <algorithm>
#include <cstdio>
#include <iostream>

using namespace mthesis;

int main()
{
    unsigned M = 30;
    unsigned N = 500;
    double oversampling = 30;

    auto sig = SynthExponentials(M, N, oversampling);

    auto parms = GPOFParams(M);
    auto ce = gpof(sig.get_y(), sig.get_d_t(), parms);

    auto y_rec = reconstruct_signal(ce, sig.get_d_t(), sig.get_N());

    std::vector<double> rel_err_db(sig.get_N());
    for (size_t n = 0; n < N; n++)
        rel_err_db[n] = calc_rel_err_db(y_rec[n], sig.get_y()[n]);

    printf("Maximum error: %.2f dB\n", *std::max_element(rel_err_db.begin(), rel_err_db.end()));

    return 0;
}
