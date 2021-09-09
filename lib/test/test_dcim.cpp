#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis.hpp>
BOOST_TEST_DONT_PRINT_LOG_VALUE( mthesis::EmMode ) // To avoid compiler error.

#include <random>

BOOST_AUTO_TEST_SUITE(DCIM)

using namespace mthesis;

EmMode mode_vals[] = {EmMode::TE, EmMode::TM};

cmplx rand_param() // Material parameters.
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_re(1, 15);
    std::uniform_real_distribution<real> dist_im(-1.0, 0.0);
    return cmplx(dist_re(gen), dist_im(gen));
}

BOOST_DATA_TEST_CASE(three_level_v1_near_range,
                     boost::unit_test::data::make(mode_vals) *
                     boost::unit_test::data::make(arma::linspace(0.1, 20, 8)) *
                     boost::unit_test::data::make(arma::linspace(0.0, 20, 8)),
                     mode, rho_by_lambda_0, z_by_lambda_0)
{
    FrequencyDomain fd(1e9);
    Medium ground(fd, rand_param(), rand_param());
    HalfSpace hs(fd, ground);

    real z_ = 0.1 * fd.lambda_0;
    VectorR3 r_ = {0, 0, z_};

    real nu = 0;
    bool direct_term = false;
    real rho = rho_by_lambda_0 * fd.lambda_0;
    real z = z_by_lambda_0 * fd.lambda_0;

    auto si = gf::scalar::layered_media::get_sommerfeld_integral(
                hs, nu, mode, direct_term);

    auto val_ref = si.eval_si_along_sip(rho, z, z_);

    auto dcim_3lv1 = dcim::ThreeLevelV2(si);
    auto val_dcim_3lv1 = dcim_3lv1.get_spatial_gf(z, z_, rho);

    auto rel_err_db = calc_rel_err_db(val_dcim_3lv1, val_ref);

    BOOST_CHECK(rel_err_db < -40);
}

BOOST_DATA_TEST_CASE(three_level_v2_near_range,
                     boost::unit_test::data::make(mode_vals) *
                     boost::unit_test::data::make(arma::linspace(0.1, 20, 8)) *
                     boost::unit_test::data::make(arma::linspace(0.0, 20, 8)),
                     mode, rho_by_lambda_0, z_by_lambda_0)
{
    FrequencyDomain fd(1e9);
    Medium ground(fd, rand_param(), rand_param());
    HalfSpace hs(fd, ground);

    real z_ = 0.1 * fd.lambda_0;
    VectorR3 r_ = {0, 0, z_};

    real nu = 0;
    bool direct_term = false;
    real rho = rho_by_lambda_0 * fd.lambda_0;
    real z = z_by_lambda_0 * fd.lambda_0;

    auto si = gf::scalar::layered_media::get_sommerfeld_integral(
                hs, nu, mode, direct_term);

    auto val_ref = si.eval_si_along_sip(rho, z, z_);

    auto dcim_3lv2 = dcim::ThreeLevelV2(si);
    auto val_dcim_3lv2 = dcim_3lv2.get_spatial_gf(z, z_, rho);

    auto rel_err_db = calc_rel_err_db(val_dcim_3lv2, val_ref);

    BOOST_CHECK(rel_err_db < -40);
}

BOOST_AUTO_TEST_SUITE_END()
