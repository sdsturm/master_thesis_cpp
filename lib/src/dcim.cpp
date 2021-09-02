#include <mthesis/detail/dcim.hpp>

#include <cassert>

namespace mthesis::dcim::threelevelv2
{
    real calc_T_03(real k_0, real k_rho_max_3)
    {
        return 1.0 / std::sqrt(1.0 - std::pow(k_rho_max_3 / k_0, 2)) - 1.0;
    }

    real calc_T_02(real k_0, real k_rho_max_2, real T_03)
    {
        real a = 1 / (1 + T_03);
        cmplx arg = std::pow(k_rho_max_2 / k_0, 2) - 1;
        auto ans = std::sqrt(arg) / a;
        assert(ans.imag() == 0.0);
        return ans.real();
    }

    real calc_T_01(real k_0, real k_rho_max_1, real T_02, real T_03)
    {
        real a = T_02 / (1 + T_03);
        return std::sqrt(std::pow(k_rho_max_1 / k_0, 2) - 1) / a;
    }

    utils::SamplingPath calc_C_3(real k_0, real T_03, int N)
    {
        auto calc_k_z = [&](real t)
        { return k_0 * (1 - t / (T_03 + 1)); };

        real d_t = utils::get_d_t(T_03, N);

        std::vector<cmplx> k_z_vals(N);
        for (int n = 0; n < N; n++)
            k_z_vals[n] = calc_k_z(n * d_t);

        return utils::SamplingPath(k_0, k_z_vals, d_t);
    }

    utils::SamplingPath calc_C_2(real k_0, real T_02, real T_03, int N)
    {
        using std::complex_literals::operator""i;

        auto calc_k_z = [&](real t)
        { return k_0 * (1 / (T_03 + 1)) * (-1.0i * t + (1 - t / T_02)); };

        auto d_t = utils::get_d_t(T_02, N);

        std::vector<cmplx> k_z_vals(N);
        for (int n = 0; n < N; n++)
            k_z_vals[n] = calc_k_z(n * d_t);

        return utils::SamplingPath(k_0, k_z_vals, d_t);
    }

    utils::SamplingPath calc_C_1(real k_0, real T_01, real T_02, real T_03, int N)
    {
        using std::complex_literals::operator""i;

        auto calc_k_z = [&](real t)
        { return -1.0i * k_0 * (T_02 / (T_03 + 1) + t); };

        auto d_t = utils::get_d_t(T_01, N);

        std::vector<cmplx> k_z_vals(N);
        for (int n = 0; n < N; n++)
            k_z_vals[n] = calc_k_z(n * d_t);

        return utils::SamplingPath(k_0, k_z_vals, d_t);
    }

    utils::ce_vec calc_coeffs_1(const utils::ce_vec &ce_in,
                                real k_0,
                                real T_02,
                                real T_03)
    {
        using std::complex_literals::operator""i;

        auto calc_alpha = [&](cmplx beta)
        { return beta / (1.0i * k_0); };

        auto calc_a = [&](cmplx b, cmplx alpha)
        { return b * std::exp(-1.0i * k_0 * alpha * T_02 / (T_03 + 1)); };

        utils::ce_vec ce_out(ce_in.size());
        for (size_t n = 0; n < ce_in.size(); n++)
        {
            auto alpha = calc_alpha(ce_in[n].back());
            auto a = calc_a(ce_in[n].front(), alpha);
            ce_out[n] = {a, alpha};
        }

        return ce_out;
    }

    utils::ce_vec calc_coeffs_2(const utils::ce_vec &ce_in,
                                real k_0,
                                real T_02,
                                real T_03)
    {
        using std::complex_literals::operator""i;

        auto calc_alpha = [&](cmplx beta)
        { return beta * T_02 * (T_03 + 1) / (k_0 * (1.0 + 1.0i * T_02)); };

        auto calc_a = [&](cmplx b, cmplx alpha)
        { return b * std::exp(k_0 * alpha / (T_03 + 1)); };

        utils::ce_vec ce_out(ce_in.size());
        for (size_t n = 0; n < ce_in.size(); n++)
        {
            auto alpha = calc_alpha(ce_in[n].back());
            auto a = calc_a(ce_in[n].front(), alpha);
            ce_out[n] = {a, alpha};
        }

        return ce_out;
    }

    utils::ce_vec calc_coeffs_3(const utils::ce_vec &ce_in,
                                real k_0,
                                real T_03)
    {
        using std::complex_literals::operator""i;

        auto calc_alpha = [&](cmplx beta)
        { return beta * (T_03 + 1) / k_0; };

        auto calc_a = [&](cmplx b, cmplx alpha)
        { return b * std::exp(k_0 * alpha); };

        utils::ce_vec ce_out(ce_in.size());
        for (size_t n = 0; n < ce_in.size(); n++)
        {
            auto alpha = calc_alpha(ce_in[n].back());
            auto a = calc_a(ce_in[n].front(), alpha);
            ce_out[n] = {a, alpha};
        }

        return ce_out;
    }

    std::vector<utils::ce_vec> three_level_v2(const si::SpectralGF &gf)
    {
        auto k_0 = utils::get_k_0(gf.lm);
        auto k_max = utils::find_k_max(gf.lm);

        // Set up sampling paths.
        real k_rho_max_3 = 0.8 * k_0; // empirical
        auto T_03 = calc_T_03(k_0, k_rho_max_3);
        int N_3 = 100;

        real k_rho_max_2 = 1.2 * k_max; // empirical
        auto T_02 = calc_T_02(k_0, k_rho_max_2, T_03);
        int N_2 = 100;

        real k_rho_max_1 = 300 * k_max; // empirical
        auto T_01 = calc_T_01(k_0, k_rho_max_1, T_02, T_03);
        int N_1 = 100;

        std::vector<utils::SamplingPath> sp;
        sp.push_back(calc_C_1(k_0, T_01, T_02, T_03, N_1));
        sp.push_back(calc_C_2(k_0, T_02, T_03, N_2));
        sp.push_back(calc_C_3(k_0, T_03, N_3));

        // Set up coefficient transforms.
        auto f1 = [=](const utils::ce_vec &ce_in)
        { return calc_coeffs_1(ce_in, k_0, T_02, T_03); };

        auto f2 = [=](const utils::ce_vec &ce_in)
        { return calc_coeffs_2(ce_in, k_0, T_02, T_03); };

        auto f3 = [=](const utils::ce_vec &ce_in)
        { return calc_coeffs_3(ce_in, k_0, T_03); };

        std::vector<utils::ct_fun> ct_funs;
        ct_funs.push_back(f1);
        ct_funs.push_back(f2);
        ct_funs.push_back(f3);

        auto ce_levels = utils::algo(gf, sp, ct_funs);

        return ce_levels;
    }

} // namespace mthesis::dcim::threelevelv2
