#include <mthesis/dcim/three_level_v2.hpp>

#include <cassert>

namespace mthesis::dcim {

ThreeLevelV2::ThreeLevelV2(const SommerfeldIntegral &si) : DCIM(si)
{
#if 0
    // Set up sampling paths.
    real k_rho_max_3 = 0.8 * k_0; // empirical
    auto T_3 = calc_T_3(k_rho_max_3);
    int N_3 = 100;

    real k_rho_max_2 = 1.2 * k_max; // empirical
    auto T_2 = calc_T_2(k_rho_max_2, T_3);
    int N_2 = 100;

    real k_rho_max_1 = 300 * k_max; // empirical
    auto T_1 = calc_T_1(k_rho_max_1, T_2, T_3);
    int N_1 = 100;
#else
    // Set up sampling paths.
    real k_rho_max_3 = 0.8 * k_0; // empirical
    auto T_3 = calc_T_3(k_rho_max_3);
    int N_3 = 200;

    real k_rho_max_2 = 1.2 * k_max; // empirical
    auto T_2 = calc_T_2(k_rho_max_2, T_3);
    int N_2 = 200;

    real k_rho_max_1 = 400 * k_max; // empirical
    auto T_1 = calc_T_1(k_rho_max_1, T_2, T_3);
    int N_1 = 300;
#endif

    sampling_paths.push_back(calc_sp_1(T_1, T_2, T_3, N_1));
    sampling_paths.push_back(calc_sp_2(T_2, T_3, N_2));
    sampling_paths.push_back(calc_sp_3(T_3, N_3));

    // Set up coefficient transforms.
    auto f1 = [=](const CeVec &ce_in)
    {
        return calc_coeffs_1(ce_in, T_2, T_3);
    };

    auto f2 = [=](const CeVec &ce_in)
    {
        return calc_coeffs_2(ce_in, T_2, T_3);
    };

    auto f3 = [=](const CeVec &ce_in)
    {
        return calc_coeffs_3(ce_in, T_3);
    };

    ct_funs.push_back(f1);
    ct_funs.push_back(f2);
    ct_funs.push_back(f3);
}

real ThreeLevelV2::calc_T_3(real k_rho_max_3) const
{
    return 1.0 / std::sqrt(1.0 - std::pow(k_rho_max_3 / k_0, 2)) - 1.0;
}

real ThreeLevelV2::calc_T_2(real k_rho_max_2, real T_3) const
{
    real a = 1 / (1 + T_3);
    cmplx arg = std::pow(k_rho_max_2 / k_0, 2) - 1;
    auto ans = std::sqrt(arg) / a;
    assert(ans.imag() == 0.0);
    return ans.real();
}

real ThreeLevelV2::calc_T_1(real k_rho_max_1, real T_2, real T_3) const
{
    real a = T_2 / (1 + T_3);
    return std::sqrt(std::pow(k_rho_max_1 / k_0, 2) - 1.0) - a;
}

SamplingPath ThreeLevelV2::calc_sp_3(real T_3, int N_3) const
{
    auto calc_k_z = [&](real t)
    {
        return k_0 * (1 - t / (T_3 + 1));
    };

    real d_t = utils::get_d_t(T_3, N_3);

    std::vector<cmplx> k_z_vals(N_3);
    for (int n = 0; n < N_3; n++) {
        k_z_vals[n] = calc_k_z(n * d_t);
    }

    return SamplingPath(k_0, k_z_vals, d_t);
}

SamplingPath ThreeLevelV2::calc_sp_2(real T_2, real T_3, int N_2) const
{
    using std::complex_literals::operator""i;

    auto calc_k_z = [&](real t)
    {
        return k_0 * (1 / (T_3 + 1)) * (-1.0i * t + (1 - t / T_2));
    };

    auto d_t = utils::get_d_t(T_2, N_2);

    std::vector<cmplx> k_z_vals(N_2);
    for (int n = 0; n < N_2; n++) {
        k_z_vals[n] = calc_k_z(n * d_t);
    }

    return SamplingPath(k_0, k_z_vals, d_t);
}

SamplingPath ThreeLevelV2::calc_sp_1(real T_1,
                                     real T_2,
                                     real T_3,
                                     int N_1) const
{
    using std::complex_literals::operator""i;

    auto calc_k_z = [&](real t)
    {
        return -1.0i * k_0 * (T_2 / (T_3 + 1.0) + t);
    };

    auto d_t = utils::get_d_t(T_1, N_1);

    std::vector<cmplx> k_z_vals(N_1);
    for (int n = 0; n < N_1; n++) {
        k_z_vals[n] = calc_k_z(n * d_t);
    }

    return SamplingPath(k_0, k_z_vals, d_t);
}

CeVec ThreeLevelV2::calc_coeffs_1(const CeVec &ce_in, real T_2, real T_3) const
{
    using std::complex_literals::operator""i;

    auto calc_alpha = [&](cmplx beta)
    {
        return beta / (1.0i * k_0);
    };

    auto calc_a = [&](cmplx b, cmplx alpha)
    {
        return b * std::exp(-1.0i * k_0 * alpha * T_2 / (T_3 + 1));
    };

    CeVec ce_out;
    for (size_t n = 0; n < ce_in.size(); n++) {
        auto alpha = calc_alpha(ce_in[n].exp);
        auto a = calc_a(ce_in[n].amp, alpha);
        ce_out.emplace_back(a, alpha);
    }

    return ce_out;
}

CeVec ThreeLevelV2::calc_coeffs_2(const CeVec &ce_in, real T_2, real T_3) const
{
    using std::complex_literals::operator""i;

    auto calc_alpha = [&](cmplx beta)
    {
        return beta * T_2 * (T_3 + 1) / (k_0 * (1.0 + 1.0i * T_2));
    };

    auto calc_a = [&](cmplx b, cmplx alpha)
    {
        return b * std::exp(k_0 * alpha / (T_3 + 1));
    };

    CeVec ce_out;
    for (size_t n = 0; n < ce_in.size(); n++) {
        auto alpha = calc_alpha(ce_in[n].exp);
        auto a = calc_a(ce_in[n].amp, alpha);
        ce_out.emplace_back(a, alpha);
    }

    return ce_out;
}

CeVec ThreeLevelV2::calc_coeffs_3(const CeVec &ce_in, real T_3) const
{
    using std::complex_literals::operator""i;

    auto calc_alpha = [&](cmplx beta)
    {
        return beta * (T_3 + 1) / k_0;
    };

    auto calc_a = [&](cmplx b, cmplx alpha)
    {
        return b * std::exp(k_0 * alpha);
    };

    CeVec ce_out;
    for (size_t n = 0; n < ce_in.size(); n++) {
        auto alpha = calc_alpha(ce_in[n].exp);
        auto a = calc_a(ce_in[n].amp, alpha);
        ce_out.emplace_back(a, alpha);
    }

    return ce_out;
}

} // namespace mthesis::dcim
