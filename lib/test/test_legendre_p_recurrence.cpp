#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis/legendre_p_recurrence.hpp>
#include <boost/math/special_functions/legendre.hpp>

BOOST_DATA_TEST_CASE(legendre_p_recurrence,
                     boost::unit_test::data::xrange<unsigned>(0, 30) *
                     boost::unit_test::data::xrange<double>(-1.0, 1.0, 0.1),
                     nu, x)
{
    auto ref = boost::math::legendre_p(nu, x);
    auto num = mthesis::legendre_p_recurrence(nu, x);
    BOOST_CHECK_CLOSE(num, ref, 1e-10);
}
