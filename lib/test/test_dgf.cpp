#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis/gf/dyadic.hpp>

#include <random>

BOOST_AUTO_TEST_SUITE(DGF)

using namespace mthesis;

BOOST_DATA_TEST_CASE(free_space_reciporcity,
                     boost::unit_test::data::xrange(5), n_run) // Multiple runs.
{
    (void)n_run; // Hush unused variable warning;

    FrequencyDomain fd(1e9);
    Medium medium = Vacuum(fd);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_pos(-100 * fd.lambda_0,
                                                  +100 * fd.lambda_0);

    VectorR3 r = {dist_pos(gen), dist_pos(gen), dist_pos(gen)};
    VectorR3 r_ = {dist_pos(gen), dist_pos(gen), dist_pos(gen)};

    DyadC3 d1, d2;
    real tol = 1e-10; // In percent.

    // Reciprocity relations: (10) - (12) in Michalski2005.
    d1 = gf::dyadic::free_space::G_EJ(medium, r, r_);
    d2 = arma::strans(gf::dyadic::free_space::G_EJ(medium, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }

    d1 = gf::dyadic::free_space::G_HM(medium, r, r_);
    d2 = arma::strans(gf::dyadic::free_space::G_HM(medium, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }

    d1 = gf::dyadic::free_space::G_EM(medium, r, r_);
    d2 = -arma::strans(gf::dyadic::free_space::G_HJ(medium, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }
}

BOOST_DATA_TEST_CASE(free_space_reference,
                     boost::unit_test::data::xrange(5), n_run) // Multiple runs.
{
    (void)n_run; // Hush unused variable warning;

    FrequencyDomain fd(1e9);
    Medium medium = Vacuum(fd);
    LayeredMedium lm(fd, {-INFINITY, INFINITY}, {medium}, BC::open, BC::open);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_pos(-10 * fd.lambda_0,
                                                  +10 * fd.lambda_0);

    VectorR3 r = {dist_pos(gen), dist_pos(gen), dist_pos(gen)};
    VectorR3 r_ = {dist_pos(gen), dist_pos(gen), dist_pos(gen)};

    DyadC3 d1, d2;
    real tol = 1e-4; // In percent.

    d1 = gf::dyadic::free_space::G_EJ(medium, r, r_);
    d2 = gf::dyadic::layered_media::G_EJ(lm, r, r_);
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }

    d1 = gf::dyadic::free_space::G_HM(medium, r, r_);
    d2 = gf::dyadic::layered_media::G_HM(lm, r, r_);
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }

    d1 = gf::dyadic::free_space::G_EM(medium, r, r_);
    d2 = gf::dyadic::layered_media::G_EM(lm, r, r_);
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }

    d1 = gf::dyadic::free_space::G_HJ(medium, r, r_);
    d2 = gf::dyadic::layered_media::G_HJ(lm, r, r_);
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }
}

BOOST_DATA_TEST_CASE(half_space_reciprocity,
                     boost::unit_test::data::xrange(5), n_run) // Multiple runs.
{
    (void)n_run; // Hush unused variable warning;
    using std::complex_literals::operator""i;

    FrequencyDomain fd(1e9);
    Medium ground(fd, 3.0 - 0.2i, 1.0);
    HalfSpace hs(fd, ground);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_z(1e-6 * fd.lambda_0,
                                                20.0 * fd.lambda_0);
    std::uniform_real_distribution<real> dist_rho(-20 * fd.lambda_0,
                                                  +20.0 * fd.lambda_0);

    // Make sure to cross the interface.
    real z = -dist_z(gen);
    real z_ = dist_z(gen);
    VectorR3 r =  {dist_rho(gen), dist_rho(gen), z};
    VectorR3 r_ = {dist_rho(gen), dist_rho(gen), z_};

    DyadC3 d1, d2;
    real tol = 1e-10; // In percent.

    // Reciprocity relations: (10) - (12) in Michalski2005.
    d1 = gf::dyadic::layered_media::G_EJ(hs, r, r_);
    d2 = arma::strans(gf::dyadic::layered_media::G_EJ(hs, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }

    d1 = gf::dyadic::layered_media::G_HM(hs, r, r_);
    d2 = arma::strans(gf::dyadic::layered_media::G_HM(hs, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }

    d1 = gf::dyadic::layered_media::G_EM(hs, r, r_);
    d2 = -arma::strans(gf::dyadic::layered_media::G_HJ(hs, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            BOOST_CHECK_CLOSE(d1(i, j).real(), d2(i, j).real(), tol);
            BOOST_CHECK_CLOSE(d1(i, j).imag(), d2(i, j).imag(), tol);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
