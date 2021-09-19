#include <mthesis.hpp>
#include <mthesis/dcim/synth_exponentials.hpp>

#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    unsigned M = 10;
    unsigned N = 201;
    double oversampling = 10;

    auto sig = dcim::gpof::SynthExponentials(M, N, oversampling);

    std::vector<cmplx> y_sub(sig.y.begin(), sig.y.begin() + N / 2);

    auto params = dcim::gpof::Params();
    params.set_M(M);
//    auto ce = dcim::gpof::gpof(sig.y, sig.d_t, params);
    auto ce = dcim::gpof::gpof(y_sub, sig.d_t, params);

    auto y_rec = dcim::gpof::reconstruct_signal(ce, sig.d_t, sig.N);

    printf("t y_true_re y_true_im y_rec_re y_rec_im rel_err_db\n");
    for (size_t n = 0; n < N; n++) {
        printf("%.6e %.6e %.6e %.6e %.6e %.6f\n",
               sig.t[n],
               sig.y[n].real(),
               sig.y[n].imag(),
               y_rec[n].real(),
               y_rec[n].imag(),
               calc_rel_err_db(y_rec[n], sig.y[n]));
    }

    return 0;
}
