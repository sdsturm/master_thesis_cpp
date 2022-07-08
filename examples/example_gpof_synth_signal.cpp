#include <mthesis.hpp>
#include <mthesis/dcim/synth_exponentials.hpp>

#include <armadillo>
#include "./../submodules/gnuplot-iostream/gnuplot-iostream.h"

#include <algorithm>
#include <iostream>
#include <cstdio>

using namespace mthesis;
using namespace mthesis::dcim;

int main()
{
    unsigned M = 10;
    unsigned N = 100;
    double oversampling = 10;

    auto sig = gpof::SynthExponentials(M, N, oversampling);

    auto params = gpof::Params();
    params.set_M(M);
    auto ce = gpof::gpof(sig.y, sig.d_t, params);

    auto y_rec = gpof::reconstruct_signal(ce, sig.d_t, sig.N);

    std::vector<double> rel_err_db(sig.N);
    std::vector<real> y_true_re(sig.N), y_true_im(sig.N), y_gpof_re(sig.N), y_gpof_im(sig.N);
    for (size_t n = 0; n < N; n++)
    {
        rel_err_db[n] = calc_rel_err_db(y_rec[n], sig.y[n]);
        y_true_re[n] = sig.y[n].real();
        y_true_im[n] = sig.y[n].imag();
        y_gpof_re[n] = y_rec[n].real();
        y_gpof_im[n] = y_rec[n].imag();
    }

    printf("Maximum error: %.2f dB\n", *std::max_element(rel_err_db.begin(), rel_err_db.end()));

    Gnuplot gp;
    gp << "set multiplot layout 3,1\n";
    gp << "set grid\n";
    gp << "plot '-' with lines title 'Re(y_true)', '-' with lines title 'Im(y_true)'\n";
    gp.send1d(boost::make_tuple(sig.t, y_true_re));
    gp.send1d(boost::make_tuple(sig.t, y_true_im));
    gp << "plot '-' with lines title 'Re(y_gpof)', '-' with lines title 'Im(y_gpof)'\n";
    gp.send1d(boost::make_tuple(sig.t, y_gpof_re));
    gp.send1d(boost::make_tuple(sig.t, y_gpof_im));
    gp << "plot '-' with lines title 'relative error in dB'\n";
    gp.send1d(boost::make_tuple(sig.t, rel_err_db));

    return 0;
}
