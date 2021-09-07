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

real rand_frequency()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist(1.0, 1e9);
    return dist(gen);
}

cmplx rand_k_rho(const FrequencyDomain &fd)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_re(0.0, 2.0 * fd.k_0);
    std::uniform_real_distribution<real> dist_im(0.0, 1.0 * fd.k_0);
    return cmplx(dist_re(gen), dist_im(gen));
}

cmplx rand_param() // Material parameters.
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_re(1, 15);
    std::uniform_real_distribution<real> dist_im(-1.0, 0.0);
    return cmplx(dist_re(gen), dist_im(gen));
}

BOOST_DATA_TEST_CASE(generic_spectral_gf_reciprocity_in_half_space_same_layer,
                     boost::unit_test::data::make(mode_vals) *
                     boost::unit_test::data::xrange(15), // Multiple runs.
                     mode, n_run)
{
    // Note: generic spectral GF is "like" V_i and I_v
    //       -> check using reciprocity relation (38) in Michalski2005.

    (void)n_run;	// Hush the unused variable warning.

    FrequencyDomain fd(rand_frequency());
    auto ground = Medium(fd, rand_param(), rand_param());
    auto lm = HalfSpace(fd, ground);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_pos(0.0, 4.0 * fd.lambda_0);

    // Note: generic spectral GF is only reciprocal for m == n.
    real z =  dist_pos(gen);
    real z_ = dist_pos(gen);

    cmplx k_rho = rand_k_rho(fd);

    bool direct_term = true;

    using gf::scalar::layered_media::generic_spectral_gf;
    cmplx val_1 = generic_spectral_gf(lm, k_rho, z, z_, mode, direct_term);
    cmplx val_2 = generic_spectral_gf(lm, k_rho, z_, z, mode, direct_term);

    real tol_in_percent = 1e-6;
    BOOST_CHECK_CLOSE(val_1.real(), val_2.real(), tol_in_percent);
    BOOST_CHECK_CLOSE(val_1.imag(), val_2.imag(), tol_in_percent);
}

BOOST_DATA_TEST_CASE(sommerfeld_identity,
                     boost::unit_test::data::make(mode_vals) *
                     boost::unit_test::data::make(arma::logspace(-1, 3, 20)) *
                     boost::unit_test::data::make(arma::logspace(-1, 3, 20)),
                     mode, x, z)
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
    VectorR3 r = {x, 0, z};
    r *= fd.lambda_0;

    cmplx ref = gf::scalar::free_space::G_0(vacuum, r, r_);

    cmplx num = si.eval_si_along_sip(r, r_);

    real rel_err_db = calc_rel_err_db(num, ref);

    BOOST_CHECK(rel_err_db < -100);
}

BOOST_AUTO_TEST_SUITE_END()
