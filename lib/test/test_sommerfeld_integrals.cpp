#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis.hpp>
BOOST_TEST_DONT_PRINT_LOG_VALUE( mthesis::EmMode ) // To avoid compiler error.

#include <armadillo>

#include <random>

BOOST_AUTO_TEST_SUITE(sommerfeld_integrals)

using namespace mthesis;

EmMode mode_vals[] = {EmMode::TE, EmMode::TM};

BOOST_DATA_TEST_CASE(sommerfeld_identity,
                     boost::unit_test::data::make(mode_vals),
                     mode)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_frequency(1.0, 1e10);

    auto fd = FrequencyDomain(dist_frequency(gen));
    auto vacuum = Vacuum(fd);
    auto lm = FreeSpace(fd);
    bool direct_term = true;
    real nu = 0;
    auto si = gf::scalar::layered_media::get_sommerfeld_integral(lm, nu, mode,
                                                                 direct_term);

    VectorR3 r_ = {0, 0, 0};
    r_ *= fd.lambda_0;

    arma::vec x_vals = arma::logspace(-1, 3, 20) * fd.lambda_0;
    arma::vec z_vals = arma::logspace(-1, 3, 20) * fd.lambda_0;

    for (auto &x : x_vals)
    {
        for (auto &z : z_vals)
        {
            VectorR3 r = {x, 0, z};
            cmplx ref = gf::scalar::free_space::G_0(vacuum, r, r_);

            // Multiple runs for timer.
            cmplx num = si.eval_si_along_sip(r, r_);

            real rel_err_db = calc_rel_err_db(num, ref);

            BOOST_CHECK(rel_err_db < -100);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
