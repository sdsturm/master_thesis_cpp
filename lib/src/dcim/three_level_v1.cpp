#include <mthesis/dcim/three_level_v1.hpp>

namespace mthesis::dcim {

ThreeLevelV1::ThreeLevelV1(const SommerfeldIntegral &si) : DCIM(si)
{
    // Set up sampling paths.
    real T_3 = 0.01;
    real T_3_ = 0.01;
    int N_3 = 100;

    real k_rho_max_2 = 1.2 * k_max;
    auto T_2 = calc_T_2(k_rho_max_2);
    int N_2 = 100;

    real k_rho_max_1 = 300 * k_max;
    auto T_1 = calc_T_1(k_rho_max_1);
    int N_1 = 100;

    sampling_paths.push_back(calc_sp_1(T_1, T_2, N_1));
    sampling_paths.push_back(calc_sp_2(T_2, T_3, T_3_, N_2));
    sampling_paths.push_back(calc_sp_3(T_3, T_3_, N_3));

    // Set up coefficient transforms.
    auto f1 = [=](const CeVec &ce_in)
    {
        return calc_coeffs_1(ce_in, T_1, T_2);
    };

    auto f2 = [=](const CeVec &ce_in)
    {
        return calc_coeffs_2(ce_in, T_2, T_3, T_3_);
    };

    auto f3 = [=](const CeVec &ce_in)
    {
        return calc_coeffs_3(ce_in, T_3, T_3_);
    };

    ct_funs.push_back(f1);
    ct_funs.push_back(f2);
    ct_funs.push_back(f3);
}

real ThreeLevelV1::calc_T_2(real k_rho_max_2) const
{
    return std::sqrt(std::pow(k_rho_max_2 / k_0, 2) - 1.0);
}

real ThreeLevelV1::calc_T_1(real k_rho_max_1) const
{
    return std::sqrt(std::pow(k_rho_max_1 / k_0, 2) - 1.0);
}

SamplingPath ThreeLevelV1::calc_sp_3(real T_3, real T_3_, int N_3) const
{
    using std::complex_literals::operator""i;

    auto calc_k_z = [&](real t)
    {
        return k_0 * (1.0 - (1.0 - T_3_ + 1.0i * T_3) / T_3 * t);
    };

    real d_t = utils::get_d_t(T_3, N_3);

    std::vector<cmplx> k_z_vals(N_3);
    for (int n = 0; n < N_3; n++) {
        k_z_vals[n] = calc_k_z(n * d_t);
    }

    return SamplingPath(k_0, k_z_vals, d_t);
}

SamplingPath ThreeLevelV1::calc_sp_2(real T_2, real T_3, real T_3_,
                                     int N_2) const
{
    using std::complex_literals::operator""i;

    auto calc_k_z = [&](real t)
    {
        return k_0 * ((T_3_ - 1.0i * T_3) -
                      (T_3_ + 1.0i * (T_2 - T_3)) / T_2 * t);
    };

    auto d_t = utils::get_d_t(T_2, N_2);

    std::vector<cmplx> k_z_vals(N_2);
    for (int n = 0; n < N_2; n++) {
        k_z_vals[n] = calc_k_z(n * d_t);
    }

    return SamplingPath(k_0, k_z_vals, d_t);
}

SamplingPath ThreeLevelV1::calc_sp_1(real T_1, real T_2, int N_1) const
{
    using std::complex_literals::operator""i;

    auto calc_k_z = [&](real t)
    {
        return -1.0i * k_0 * (T_2 + (T_1 - T_2) / T_1 * t);
    };

    auto d_t = utils::get_d_t(T_1, N_1);

    std::vector<cmplx> k_z_vals(N_1);
    for (int n = 0; n < N_1; n++) {
        k_z_vals[n] = calc_k_z(n * d_t);
    }

    return SamplingPath(k_0, k_z_vals, d_t);
}

CeVec ThreeLevelV1::calc_coeffs_1(const CeVec &ce_in, real T_1, real T_2) const
{
    using std::complex_literals::operator""i;

    auto calc_alpha = [&](cmplx beta)
    {
        return beta * T_1 / (1.0i * k_0 * (T_1 - T_2));
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

CeVec ThreeLevelV1::calc_coeffs_2(const CeVec &ce_in, real T_2, real T_3,
                                  real T_3_) const
{
    using std::complex_literals::operator""i;

    auto calc_alpha = [&](cmplx beta)
    {
        return beta * T_2 / (k_0 * (1.0i * (T_2 - T_3) + T_3_));
    };

    auto calc_a = [&](cmplx b, cmplx alpha)
    {
        return  b * std::exp(k_0 * alpha * (T_3_ - 1.0i * T_3));
    };

    CeVec ce_out;
    for (size_t n = 0; n < ce_in.size(); n++) {
        auto alpha = calc_alpha(ce_in[n].exp);
        auto a = calc_a(ce_in[n].amp, alpha);
        ce_out.emplace_back(a, alpha);
    }

    return ce_out;
}

CeVec ThreeLevelV1::calc_coeffs_3(const CeVec &ce_in, real T_3, real T_3_) const
{
    using std::complex_literals::operator""i;

    auto calc_alpha = [&](cmplx beta)
    {
        return beta * T_3 / (k_0 * (1.0 - T_3_ + 1.0i * T_3));
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
