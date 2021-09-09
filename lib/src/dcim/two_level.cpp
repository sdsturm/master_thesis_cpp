#include <mthesis/dcim/two_level.hpp>

namespace mthesis::dcim {

TwoLevel::TwoLevel(const SommerfeldIntegral &si) : DCIM(si)
{
    // Set up sampling paths.
    real k_rho_max_2 = 1.5 * k_max;
    real T_2 = std::sqrt(std::pow(k_rho_max_2 / real(k_0), 2) - 1.0);
    int N_2 = 100;

    real k_rho_max_1 = 300 * k_max;
    real T_1 = std::sqrt(std::pow(k_rho_max_1 / real(k_0), 2) - 1.0) - T_2;
    int N_1 = 100;

    sampling_paths.push_back(calc_sp_1(T_1, T_2, N_1));
    sampling_paths.push_back(calc_sp_2(T_2, N_2));

    // Set up coefficient transforms.
    auto f1 = [=](const CeVec &ce_in)
    {
        return calc_coeffs_1(ce_in, T_2);
    };

    auto f2 = [=](const CeVec &ce_in)
    {
        return calc_coeffs_2(ce_in, T_2);
    };

    ct_funs.push_back(f1);
    ct_funs.push_back(f2);
}

SamplingPath TwoLevel::calc_sp_2(real T_2, int N_2) const
{
    using std::complex_literals::operator""i;

    auto calc_k_z = [&](real t)
    {
        return k_0 * (-1.0i * t + (1.0 - t / T_2));
    };

    auto d_t = utils::get_d_t(T_2, N_2);

    std::vector<cmplx> k_z_vals(N_2);
    for (int n = 0; n < N_2; n++) {
        k_z_vals[n] = calc_k_z(n * d_t);
    }

    return SamplingPath(k_0, k_z_vals, d_t);
}

SamplingPath TwoLevel::calc_sp_1(real T_1, real T_2, int N_1) const
{
    using std::complex_literals::operator""i;

    auto calc_k_z = [&](real t)
    {
        return -1.0 * k_0 * (T_2 + t);
    };

    auto d_t = utils::get_d_t(T_1, N_1);

    std::vector<cmplx> k_z_vals(N_1);
    for (int n = 0; n < N_1; n++) {
        k_z_vals[n] = calc_k_z(n * d_t);
    }

    return SamplingPath(k_0, k_z_vals, d_t);
}

CeVec TwoLevel::calc_coeffs_1(const CeVec &ce_in, real T_2) const
{
    using std::complex_literals::operator""i;

    auto calc_alpha = [&](cmplx beta)
    {
        return beta / (1.0i * k_0);
    };

    auto calc_a = [&](cmplx b, cmplx alpha)
    {
        return b * std::exp(-1.0i * k_0 * alpha * T_2);
    };

    CeVec ce_out;
    for (size_t n = 0; n < ce_in.size(); n++) {
        auto alpha = calc_alpha(ce_in[n].exp);
        auto a = calc_a(ce_in[n].amp, alpha);
        ce_out.emplace_back(a, alpha);
    }

    return ce_out;
}

CeVec TwoLevel::calc_coeffs_2(const CeVec &ce_in, real T_2) const
{
    using std::complex_literals::operator""i;

    auto calc_alpha = [&](cmplx beta)
    {
        return beta * T_2 / (1.0 + 1.0i * T_2) / k_0;
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
