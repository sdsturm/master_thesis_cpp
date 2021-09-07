#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis.hpp>
BOOST_TEST_DONT_PRINT_LOG_VALUE( mthesis::EmMode ) // To avoid compiler error.

#include <armadillo>

#include <random>

BOOST_AUTO_TEST_SUITE(SGF)

using namespace mthesis;

EmMode mode_vals[] = {EmMode::TE, EmMode::TM};

#if 0
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
#endif

#if 0
BOOST_DATA_TEST_CASE(generic_spectral_half_space_reciprocity,
                     boost::unit_test::data::make(mode_vals) *
                     boost::unit_test::data::xrange(5), // Multiple runs.
                     mode, n_run)
{
    // Note: generic spectral GF is "like" V_i and I_v
    //       -> check using reciprocity relation (38) in Michalski2005.

    using std::complex_literals::operator""i;
    (void)n_run; // Hush unused variable warning.

    FrequencyDomain fd(1e9);
    Medium ground(fd, 3.0 - 0.02i, 1.0);
//    Medium ground(fd, 1.0, 3.0 - 0.02i);
    HalfSpace hs(fd, ground);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_z(1e-6 * fd.lambda_0,
                                                20.0 * fd.lambda_0);
    std::uniform_real_distribution<real> dist_k_rho_re(0.0, 3.0 * fd.k_0);
    std::uniform_real_distribution<real> dist_k_rho_im(0.0, 0.2 * fd.k_0);

    // Make sure to cross the interface.
    real z = -dist_z(gen);
    real z_ = dist_z(gen);

    cmplx k_rho(dist_k_rho_re(gen), dist_k_rho_im(gen));

    bool direct_term = true;

    using gf::scalar::layered_media::generic_spectral;
    cmplx val1 = generic_spectral(hs, z, z_, k_rho, mode, direct_term);
    cmplx val2 = generic_spectral(hs, z_, z, k_rho, mode, direct_term);

    real tol = 1e-6;	// In percent.
    BOOST_CHECK_CLOSE(val1.real(), val2.real(), tol);
    BOOST_CHECK_CLOSE(val1.imag(), val2.imag(), tol);
}
#endif

BOOST_AUTO_TEST_SUITE_END()
