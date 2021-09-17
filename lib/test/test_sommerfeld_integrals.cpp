#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis/si/partition_extrapolation.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <armadillo>

#include <random>

BOOST_AUTO_TEST_SUITE(sommerfeld_integrals)

using namespace mthesis;

BOOST_DATA_TEST_CASE(bessel_zero_finder,
                     boost::unit_test::data::xrange<real>(-3, 3, 0.2) *
                     boost::unit_test::data::xrange<int>(0, 1, 1) *
                     boost::unit_test::data::xrange<int>(5, 1000, 5),
                     log_rho, nu_int, m)
{
    real rho = std::pow(10, log_rho);
    real nu = static_cast<real>(nu_int);

    real j_m = boost::math::cyl_bessel_j_zero(nu, m);
    real j_next = boost::math::cyl_bessel_j_zero(nu, m + 1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist(j_m, j_next);

    real a = dist(gen) / rho;

    real j_found = si::pe::utils::get_first_zero(nu, a, rho);
    real check_val = j_next / rho;

    BOOST_CHECK_CLOSE(j_found, check_val, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
