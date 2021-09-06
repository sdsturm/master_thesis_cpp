#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <mthesis/gf/dyadic.hpp>

#include <random>

BOOST_AUTO_TEST_SUITE(dgf)

using namespace mthesis;

BOOST_AUTO_TEST_CASE(free_space_reciporcity)
{
    FrequencyDomain fd(1e9);
    Medium medium = Vacuum(fd);

    VectorR3 r_ = {0, 0, 0};
    VectorR3 r = {1, 2, 3};
    r *= fd.lambda_0;
    r_ *= fd.lambda_0;

    DyadC3 d1, d2;

    d1 = gf::dyadic::free_space::G_EJ(medium, r, r_);
    d2 = arma::strans(gf::dyadic::free_space::G_EJ(medium, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }

    d1 = gf::dyadic::free_space::G_HM(medium, r, r_);
    d2 = arma::strans(gf::dyadic::free_space::G_HM(medium, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }

    d1 = gf::dyadic::free_space::G_EM(medium, r, r_);
    d2 = -arma::strans(gf::dyadic::free_space::G_HJ(medium, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }
}

BOOST_AUTO_TEST_CASE(free_space_reference)
{
    FrequencyDomain fd(1e9);
    Medium medium = Vacuum(fd);
    LayeredMedium lm(fd, {-INFINITY, INFINITY}, {medium}, BC::open, BC::open);

    VectorR3 r_ = {0, 0, 0};
    VectorR3 r = {1, 2, 3};
    r *= fd.lambda_0;
    r_ *= fd.lambda_0;

    DyadC3 d1, d2;

    d1 = gf::dyadic::free_space::G_EJ(medium, r, r_);
    d2 = gf::dyadic::layered_media::G_EJ(lm, r, r_);
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }

    d1 = gf::dyadic::free_space::G_EM(medium, r, r_);
    d2 = gf::dyadic::layered_media::G_EM(lm, r, r_);
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }

    d1 = gf::dyadic::free_space::G_HJ(medium, r, r_);
    d2 = gf::dyadic::layered_media::G_HJ(lm, r, r_);
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }

    d1 = gf::dyadic::free_space::G_HM(medium, r, r_);
    d2 = gf::dyadic::layered_media::G_HM(lm, r, r_);
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }
}

BOOST_AUTO_TEST_CASE(half_space_reciprocity)
{
    using std::complex_literals::operator""i;

    FrequencyDomain fd(1e9);
    Medium ground(fd, 3.0 - 0.2i, 1.0);
    HalfSpace hs(fd, ground);

    VectorR3 r_ = {0, 0, -1};
    VectorR3 r = {1, 2, 3};
    r *= fd.lambda_0;
    r_ *= fd.lambda_0;

    DyadC3 d1, d2;

    d1 = gf::dyadic::layered_media::G_EJ(hs, r, r_);
    d2 = arma::strans(gf::dyadic::layered_media::G_EJ(hs, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }

    d1 = gf::dyadic::layered_media::G_HM(hs, r, r_);
    d2 = arma::strans(gf::dyadic::layered_media::G_HM(hs, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }

    d1 = gf::dyadic::layered_media::G_EM(hs, r, r_);
    d2 = -arma::strans(gf::dyadic::layered_media::G_HJ(hs, r_, r));
    for (arma::uword i = 0; i < 3; i++) {
        for (arma::uword j = 0; j < 3; j++) {
            real rel_err_db = calc_rel_err_db(d2(i, j), d1(i, j));
            BOOST_CHECK_LE(rel_err_db, -100);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
