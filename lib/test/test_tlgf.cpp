#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis/definitions.hpp>
BOOST_TEST_DONT_PRINT_LOG_VALUE( mthesis::EmMode ) // To avoid compiler error.
#include <mthesis/solution_domain.hpp>
#include <mthesis/tlgf.hpp>

#include <armadillo>
#include <gsl/gsl_const_mksa.h>

#include <random>

BOOST_AUTO_TEST_SUITE(TLGF)

using namespace mthesis;

constexpr unsigned n_runs = 15;

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
    std::uniform_real_distribution<real> dist_im(0.0, 2.0 * fd.k_0);
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

void reciprocity_check(const LayeredMedium &lm,
                       real z,
                       real z_,
                       cmplx k_rho,
                       EmMode type)
{
    cmplx val_1, val_2;
    using namespace mthesis::tlgf;

    real tol_in_percent = 1e-6;

    auto check = [&](cmplx val_1, cmplx val_2)
    {
        BOOST_CHECK_CLOSE(val_1.real(), val_2.real(), tol_in_percent);
        BOOST_CHECK_CLOSE(val_1.imag(), val_2.imag(), tol_in_percent);
    };

    val_1 = V_i(lm, k_rho, z, z_, type);
    val_2 = V_i(lm, k_rho, z_, z, type);
    check(val_1, val_2);

    val_1 = I_v(lm, k_rho, z, z_, type);
    val_2 = I_v(lm, k_rho, z_, z, type);
    check(val_1, val_2);

    val_1 = V_v(lm, k_rho, z, z_, type);
    val_2 = -I_i(lm, k_rho, z_, z, type);
    check(val_1, val_2);

    val_1 = I_i(lm, k_rho, z, z_, type);
    val_2 = -V_v(lm, k_rho, z_, z, type);
    check(val_1, val_2);
}

BOOST_DATA_TEST_CASE(reciprocity_in_free_space,
                     boost::unit_test::data::make(mode_vals) *
                     boost::unit_test::data::xrange(n_runs), // Multiple runs.
                     type, n_run)
{
    (void)n_run;	// Hush the unused variable warning.

    FrequencyDomain fd(rand_frequency());
    LayeredMedium lm = FreeSpace(fd);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_pos(-4.0 * fd.lambda_0,
                                                  +4.0 * fd.lambda_0);

    real z = dist_pos(gen);
    real z_ = dist_pos(gen);
    cmplx k_rho = rand_k_rho(fd);

    reciprocity_check(lm, z, z_, k_rho, type);
}

BOOST_DATA_TEST_CASE(reciprocity_in_half_space,
                     boost::unit_test::data::make(mode_vals) *
                     boost::unit_test::data::xrange(n_runs), // Multiple runs.
                     type, n_run)
{
    (void)n_run;	// Hush the unused variable warning.

    FrequencyDomain fd(rand_frequency());
    auto ground = Medium(fd, rand_param(), rand_param());
    auto lm = HalfSpace(fd, ground);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_pos(1e-6 * fd.lambda_0,
                                                  4.0 * fd.lambda_0);

    // Points on opposite sides of interface.
    real z = dist_pos(gen);
    real z_ = -dist_pos(gen);

    cmplx k_rho = rand_k_rho(fd);

    reciprocity_check(lm, z, z_, k_rho, type);
}

BOOST_DATA_TEST_CASE(reciprocity_in_multilayer,
                     boost::unit_test::data::make(mode_vals) *
                     boost::unit_test::data::xrange(n_runs), // Multiple runs.
                     type, n_run)
{
    (void)n_run;	// Hush the unused variable warning.

    FrequencyDomain fd(rand_frequency());
    auto medium = Medium(fd, rand_param(), rand_param());
    std::vector<real> z_interfaces = {-INFINITY,
                                      -fd.lambda_0,
                                      fd.lambda_0,
                                      INFINITY};
    std::vector<Medium> media = {Vacuum(fd), medium, Vacuum(fd)};
    auto lm = LayeredMedium(fd, z_interfaces, media, BC::open, BC::open);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<real> dist_pos(1.1 * fd.lambda_0,
                                                  4.0 * fd.lambda_0);

    // Transmission through both interfaces.
    real z = dist_pos(gen);
    real z_ = -dist_pos(gen);

    cmplx k_rho = rand_k_rho(fd);

    reciprocity_check(lm, z, z_, k_rho, type);
}

BOOST_DATA_TEST_CASE(impedance_current_source,
                     boost::unit_test::data::make(mode_vals), mode)
{
    FrequencyDomain fd(rand_frequency());
    auto medium = Medium(fd, rand_param(), rand_param());
    std::vector<real> z_interfaces = {-INFINITY,
                                      -fd.lambda_0,
                                      fd.lambda_0,
                                      INFINITY};
    std::vector<Medium> media = {Vacuum(fd), medium, Vacuum(fd)};
    auto lm = LayeredMedium(fd, z_interfaces, media, BC::open, BC::open);

    // Soruce inside the sheet.
    real z_ = 0;

    arma::vec z_vals = arma::linspace(-4, 4, 300) * fd.lambda_0;
    cmplx k_rho = 0;

    const real Z_0 = std::sqrt(GSL_CONST_MKSA_VACUUM_PERMEABILITY /
                               GSL_CONST_MKSA_VACUUM_PERMITTIVITY);

    using namespace mthesis::tlgf;
    for (const auto &z : z_vals) {
        cmplx V = V_i(lm, k_rho, z, z_, mode);
        cmplx I = I_i(lm, k_rho, z, z_, mode);
        cmplx Z = V / I;
        if (z < *(lm.z.begin() + 1)) {
            BOOST_CHECK_CLOSE(Z.real(), -Z_0, 1e-4);
            BOOST_CHECK(std::abs(Z.imag()) < 1e-6);
        }
        if (z > *(lm.z.end() - 2)) {
            BOOST_CHECK_CLOSE(Z.real(), Z_0, 1e-4);
            BOOST_CHECK(std::abs(Z.imag()) < 1e-6);
        }
    }
}

BOOST_DATA_TEST_CASE(impedance_voltage_source,
                     boost::unit_test::data::make(mode_vals), mode)
{
    FrequencyDomain fd(rand_frequency());
    auto medium = Medium(fd, rand_param(), rand_param());
    std::vector<real> z_interfaces = {-INFINITY,
                                      -fd.lambda_0,
                                      fd.lambda_0,
                                      INFINITY};
    std::vector<Medium> media = {Vacuum(fd), medium, Vacuum(fd)};
    auto lm = LayeredMedium(fd, z_interfaces, media, BC::open, BC::open);

    // Soruce inside the sheet.
    real z_ = 0;
    arma::vec z_vals = arma::linspace(-4, 4, 300) * fd.lambda_0;
    cmplx k_rho = 0;

    const real Z_0 = std::sqrt(GSL_CONST_MKSA_VACUUM_PERMEABILITY /
                               GSL_CONST_MKSA_VACUUM_PERMITTIVITY);

    using namespace mthesis::tlgf;
    for (const auto &z : z_vals) {
        cmplx V = V_v(lm, k_rho, z, z_, mode);
        cmplx I = I_v(lm, k_rho, z, z_, mode);
        cmplx Z = V / I;
        if (z < *(lm.z.begin() + 1)) {
            BOOST_CHECK_CLOSE(Z.real(), -Z_0, 1e-4);
            BOOST_CHECK(std::abs(Z.imag()) < 1e-6);
        }
        if (z > *(lm.z.end() - 2)) {
            BOOST_CHECK_CLOSE(Z.real(), Z_0, 1e-4);
            BOOST_CHECK(std::abs(Z.imag()) < 1e-6);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
