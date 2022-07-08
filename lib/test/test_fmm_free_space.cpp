#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <mthesis/definitions.hpp>
#include <mthesis/fmm.hpp>
#include <mthesis/mom.hpp>

#include <random>

BOOST_AUTO_TEST_SUITE(FMM_free_space)

using namespace mthesis;

BOOST_DATA_TEST_CASE(plane_wave_integral,
                     boost::unit_test::data::xrange<int>(10, 50, 10), L)
{
    fmm::EwaldSphere es(L);
    std::vector<cmplx> f(es.k_hat.size(), 1.0);
    cmplx num = es.integrate(f);
    cmplx ref(4.0 * M_PI, 0.0);

    BOOST_CHECK_CLOSE(num.real(), ref.real(), 1e-8);
    BOOST_CHECK_CLOSE(num.imag(), ref.imag(), 1e-8);
}

BOOST_DATA_TEST_CASE(point_generation_and_group_identification,
                     boost::unit_test::data::xrange(7), // Run test 7 times.
                     run_numer)
{
    (void)run_numer; // Just to hush unused variable warning.

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist_ind(-100, 100);
    std::uniform_real_distribution<real> dist_w_factor(1.0, 10.0);

    FrequencyDomain fd(1e9);
    real w = fd.lambda_0 * dist_w_factor(gen);
    fmm::Params params(fd, w);

    fmm::multiindex mi = {dist_ind(gen), dist_ind(gen), dist_ind(gen)};

    unsigned N_pts = 10;
    auto pts = fmm::rand_pts_in_group(params, mi, N_pts);

    for (const auto &r : pts) {
        auto mi_found = fmm::identify_group(params, r);
        BOOST_CHECK(arma::all(mi == mi_found));
    }
}

BOOST_DATA_TEST_CASE(group_building,
                     boost::unit_test::data::xrange(7), // Run test 7 times.
                     run_number)
{
    (void)run_number; // Just to hush unused variable warning.

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist_ind(-100, 100);
    std::uniform_real_distribution<real> dist_w_factor(1.0, 10.0);

    FrequencyDomain fd(1e9);
    real w = fd.lambda_0 * dist_w_factor(gen);
    fmm::Params params(fd, w);

    fmm::multiindex mi_1 = {dist_ind(gen), dist_ind(gen), dist_ind(gen)};
    fmm::multiindex mi_2 = {dist_ind(gen), dist_ind(gen), dist_ind(gen)};
    // Make sure multiindices are different.
    while (arma::all(mi_1 == mi_2)) {
        fmm::multiindex mi_1 = {dist_ind(gen), dist_ind(gen), dist_ind(gen)};
        fmm::multiindex mi_2 = {dist_ind(gen), dist_ind(gen), dist_ind(gen)};
    }

    unsigned N_pts = 10;
    auto pts = fmm::rand_pts_in_group(params, mi_1, N_pts);
    fmm::append_pts(pts, fmm::rand_pts_in_group(params, mi_2, N_pts));

    auto groups = fmm::build_groups(params, pts);

    BOOST_CHECK(groups.size() == 2);
    BOOST_CHECK(arma::all(groups.front().mi == mi_1));
    BOOST_CHECK(arma::all(groups.back().mi == mi_2));
}

BOOST_AUTO_TEST_CASE(FMM_with_random_points)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist_ind(-100, 100);
    std::uniform_real_distribution<real> dist_frequency(1.0, 1e10);
    std::uniform_real_distribution<real> dist_w_factor(1.0, 1.5);
    std::uniform_real_distribution<real> dist_amplitide(-1.0, 1-0);

    auto rand_mi = [&]()
    {
        return fmm::multiindex {dist_ind(gen), dist_ind(gen), dist_ind(gen)};
    };

    FrequencyDomain fd(dist_frequency(gen));
    real w = fd.lambda_0 * dist_w_factor(gen);
    fmm::Params params(fd, w);

label:
    unsigned N = 30;
    std::vector<VectorR3> src_pts;
    std::vector<VectorR3> obs_pts;
    try {
        fmm::append_pts(src_pts, fmm::rand_pts_in_group(params, rand_mi(), N));
        fmm::append_pts(src_pts, fmm::rand_pts_in_group(params, rand_mi(), N));
        fmm::append_pts(src_pts, fmm::rand_pts_in_group(params, rand_mi(), N));

        fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, rand_mi(), N));
        fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, rand_mi(), N));
        fmm::append_pts(obs_pts, fmm::rand_pts_in_group(params, rand_mi(), N));

        fmm::freespace::FMM(params, src_pts, obs_pts); // To check if valid.
    }  catch (fmm::GroupSeparationError &e) {
        goto label;
    }

    fmm::freespace::FMM my_fmm(params, src_pts, obs_pts);
    mom::MoM my_mom(fd, src_pts, obs_pts);

    std::vector<cmplx> I(src_pts.size());
    for (auto &elem : I) {
        elem = cmplx(dist_amplitide(gen), dist_amplitide(gen));
    }

    auto V_fmm = my_fmm.calc_product(I);
    auto V_mom = my_mom.calc_product(I);

    real rel_err_db_max = -std::numeric_limits<real>::infinity();
    for (size_t n = 0; n < obs_pts.size(); n++) {
        real rel_err_db = calc_rel_err_db(V_fmm[n], V_mom[n]);
        if (rel_err_db > rel_err_db_max) {
            rel_err_db_max = rel_err_db;
        }
    }

    BOOST_CHECK(rel_err_db_max < -150);
}

BOOST_AUTO_TEST_SUITE_END()
