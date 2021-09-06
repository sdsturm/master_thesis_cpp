#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis/gpof.hpp>
#include <mthesis/misc/synth_exponentials.hpp>

BOOST_AUTO_TEST_SUITE(GPOF)

using namespace mthesis;

constexpr unsigned M_min = 3;
constexpr unsigned M_max = 21;
constexpr unsigned M_step = 3;

constexpr unsigned N_min = 50;
constexpr unsigned N_max = 300;
constexpr unsigned N_step = 50;

constexpr real oversampling_vals[] = {2.0001, 10, 20, 30, 40, 50};

constexpr real tol_vals[] = {1e-6, 1e-9, 1e-12};

real get_rel_err_db_max(const std::vector<cmplx> &y_reconstructed,
                        const std::vector<cmplx> &y_true)
{
    assert(y_reconstructed.size() == y_true.size());

    real rel_err_db_max = -std::numeric_limits<real>::infinity();
    for (size_t n = 0; n < y_true.size(); n++) {
        real rel_err_db = calc_rel_err_db(y_reconstructed[n], y_true[n]);
        if (rel_err_db > rel_err_db_max) {
            rel_err_db_max = rel_err_db;
        }
    }
    return rel_err_db_max;
}

BOOST_DATA_TEST_CASE(given_model_order,
                     boost::unit_test::data::xrange<unsigned>(M_min, M_max, M_step) *
                     boost::unit_test::data::xrange<unsigned>(N_min, N_max, N_step) *
                     boost::unit_test::data::make(oversampling_vals),
                     M, N, oversampling)
{
    gpof::SynthExponentials sig(M, N, oversampling);

    gpof::Params params;
    params.set_M(M);
    auto ce = gpof::gpof(sig.y, sig.d_t, params);

    auto y_rec = gpof::reconstruct_signal(ce, sig.d_t, sig.N);

    auto rel_err_db_max = get_rel_err_db_max(y_rec, sig.y);

    real check_value = -150;

    BOOST_TEST(rel_err_db_max < check_value);
}

BOOST_DATA_TEST_CASE(given_tolerance,
                     boost::unit_test::data::xrange<unsigned>(M_min, M_max, M_step) *
                     boost::unit_test::data::xrange<unsigned>(N_min, N_max, N_step) *
                     boost::unit_test::data::make(oversampling_vals) *
                     boost::unit_test::data::make(tol_vals),
                     M, N, oversampling, tol)
{
    gpof::SynthExponentials sig(M, N, oversampling);

    gpof::Params params;
    params.set_tol(tol);
    auto ce = gpof::gpof(sig.y, sig.d_t, params);

    auto y_rec = gpof::reconstruct_signal(ce, sig.d_t, sig.N);

    auto rel_err_db_max = get_rel_err_db_max(y_rec, sig.y);

    real check_value = 20.0 * std::log10(tol);
    check_value += 0.4 * std::abs(check_value); // Factor set empirically.

    BOOST_TEST(rel_err_db_max < check_value);
}

BOOST_AUTO_TEST_SUITE_END()
