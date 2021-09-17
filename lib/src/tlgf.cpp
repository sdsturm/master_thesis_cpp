#include <mthesis/tlgf.hpp>

namespace mthesis::tlgf {

cmplx V_i(const LayeredMedium &lm,
          cmplx k_rho,
          real z,
          real z_,
          EmMode type,
          bool direct_term,
          RiemannSheet sheet)
{
    auto m = lm.identify_layer(z);
    auto n = lm.identify_layer(z_);
    bool dual_solution = false;
    auto d = utils::Internals(lm, k_rho, type, dual_solution, sheet);
    return utils::V_i_base(lm, z, z_, m, n, d, direct_term);
}

cmplx I_i(const LayeredMedium &lm,
          cmplx k_rho,
          real z,
          real z_,
          EmMode type,
          bool direct_term,
          RiemannSheet sheet)
{
    auto m = lm.identify_layer(z);
    auto n = lm.identify_layer(z_);
    bool dual_solution = false;
    auto d = utils::Internals(lm, k_rho, type, dual_solution, sheet);
    return utils::I_i_base(lm, z, z_, m, n, d, direct_term);
}

cmplx I_v(const LayeredMedium &lm,
          cmplx k_rho,
          real z,
          real z_,
          EmMode type,
          bool direct_term,
          RiemannSheet sheet)
{
    auto m = lm.identify_layer(z);
    auto n = lm.identify_layer(z_);
    bool dual_solution = true;
    auto d = utils::Internals(lm, k_rho, type, dual_solution, sheet);
    return utils::V_i_base(lm, z, z_, m, n, d, direct_term);
}

cmplx V_v(const LayeredMedium &lm,
          cmplx k_rho,
          real z,
          real z_,
          EmMode type,
          bool direct_term,
          RiemannSheet sheet)
{
    auto m = lm.identify_layer(z);
    auto n = lm.identify_layer(z_);
    bool dual_solution = true;
    auto d = utils::Internals(lm, k_rho, type, dual_solution, sheet);
    return utils::I_i_base(lm, z, z_, m, n, d, direct_term);
}

cmplx generic_sgf(const LayeredMedium &lm,
                  cmplx k_rho,
                  real z,
                  real z_,
                  EmMode type,
                  bool direct_term,
                  RiemannSheet sheet)
{
    using std::complex_literals::operator""i;
    auto m = lm.identify_layer(z);
    auto n = lm.identify_layer(z_);
    bool dual_solution = false;
    auto d = utils::Internals(lm, k_rho, type, dual_solution, sheet);

    // Note: this function is not reciprocal if z and z_ lie in layers with
    // 	     different material paramters.
    //       For the case m == n (z and z_ in the same layer) this function
    //       fulfills reciprocity.
    return utils::V_i_base_wo_prefac(lm, z, z_, m, n, d, direct_term) /
            (1.0i * d.k_z[n]);
}

// Hide implementation details in nested namespace utils.
namespace utils {

void calc_k_z_select_sheet(std::vector<cmplx> &k_z,
                           RiemannSheet sheet)
{
    switch (sheet) {
    case RiemannSheet::I:
        // Top layer.
        if (k_z.back().imag() > 0.0)
            k_z.back() *= -1.0;
        // Bottom layer.
        if (k_z.front().imag() > 0.0)
            k_z.front() *= -1.0;
        break;
    case RiemannSheet::II:
        // Top layer.
        if (k_z.back().imag() < 0.0)
            k_z.back() *= -1.0;
        // Bottom layer.
        if (k_z.front().imag() > 0.0)
            k_z.front() *= -1.0;
        break;
    case RiemannSheet::III:
        // Top layer.
        if (k_z.back().imag() > 0.0)
            k_z.back() *= -1.0;
        // Bottom layer.
        if (k_z.front().imag() < 0.0)
            k_z.front() *= -1.0;
        break;
    case RiemannSheet::IV:
        // Top layer.
        if (k_z.back().imag() < 0.0)
            k_z.back() *= -1.0;
        // Bottom layer.
        if (k_z.front().imag() < 0.0)
            k_z.front() *= -1.0;
        break;
    }
}

std::vector<cmplx> calc_k_z(const LayeredMedium &lm,
                            cmplx k_rho,
                            RiemannSheet sheet)
{
    std::vector<cmplx> k_z(lm.media.size());
    cmplx k_rho_square = std::pow(k_rho, 2.0);

    for (size_t n = 0; n < lm.media.size(); n++) {
        k_z[n] = std::sqrt(std::pow(lm.media[n].k, 2.0) - k_rho_square);
    }

    calc_k_z_select_sheet(k_z, sheet);

    return k_z;
}

std::vector<cmplx> calc_Z(const LayeredMedium &lm,
                          const std::vector<cmplx> &k_z,
                          EmMode type)
{
    std::vector<cmplx> Z(lm.media.size());

    switch (type) {
    case EmMode::TM:
        for (size_t n = 0; n < lm.media.size(); n++) {
            Z[n] = k_z[n] / (lm.fd.omega * lm.media[n].eps);
        }
        break;
    case EmMode::TE:
        for (size_t n = 0; n < lm.media.size(); n++) {
            Z[n] = (lm.fd.omega * lm.media[n].mu) / k_z[n];
        }
        break;
    }

    return Z;
}

std::vector<cmplx> calc_Y(const LayeredMedium &lm,
                          const std::vector<cmplx> &k_z,
                          EmMode type)
{
    auto Y = calc_Z(lm, k_z, type);
    for (auto &elem : Y) {
        elem = 1.0 / elem;
    }

    return Y;
}

std::vector<cmplx> calc_theta(const LayeredMedium &lm,
                              const std::vector<cmplx> &k_z)
{
    std::vector<cmplx> theta(lm.media.size());

    for (size_t n = 0; n < lm.media.size(); n++) {
        theta[n] = k_z[n] * lm.d[n];
    }

    return theta;
}

cmplx calc_Gamma_fresnel(const std::vector<cmplx> &Z, int i, int j)
{
    return (Z[i] - Z[j]) / (Z[i] + Z[j]);
}

std::vector<cmplx> calc_Gamma_u(const LayeredMedium &lm,
                                const std::vector<cmplx> &Z,
                                const std::vector<cmplx> &theta)
{
    std::vector<cmplx> Gamma_u(lm.media.size());

    switch (lm.top_bc) {
    case BC::open:
        Gamma_u.back() = 0.0;
        break;
    case BC::PEC:
        Gamma_u.back() = -1.0;
        break;
    case BC::PMC:
        Gamma_u.back() = 1.0;
        break;
    }

    using std::complex_literals::operator""i;
    for (int n = lm.media.size() - 2; n >= 0; n--) {
        auto t1 = calc_Gamma_fresnel(Z, n + 1, n);
        if (std::isfinite(lm.d[n + 1])) {
            cmplx t2 = Gamma_u[n + 1] * std::exp(-2.0i * theta[n + 1]);
            Gamma_u[n] = (t1 + t2) / (1.0 + t1 * t2);
        } else {
            Gamma_u[n] = t1;
        }
    }

    return Gamma_u;
}

std::vector<cmplx> calc_Gamma_d(const LayeredMedium &lm,
                                const std::vector<cmplx> &Z,
                                const std::vector<cmplx> &theta)
{
    std::vector<cmplx> Gamma_d(lm.media.size());

    switch (lm.bottom_bc) {
    case BC::open:
        Gamma_d.front() = 0.0;
        break;
    case BC::PEC:
        Gamma_d.front() = -1.0;
        break;
    case BC::PMC:
        Gamma_d.front() = 1.0;
        break;
    }

    using std::complex_literals::operator""i;
    for (size_t n = 1; n < lm.media.size(); n++) {
        auto t1 = calc_Gamma_fresnel(Z, n - 1, n);
        if (std::isfinite(lm.d[n - 1])) {
            cmplx t2 = Gamma_d[n - 1] * std::exp(-2.0i * theta[n - 1]);
            Gamma_d[n] = (t1 + t2) / (1.0 + t1 * t2);
        } else {
            Gamma_d[n] = t1;
        }
    }

    return Gamma_d;
}

Internals::Internals(const LayeredMedium &lm,
                     cmplx k_rho,
                     EmMode mode,
                     bool dual_solution,
                     RiemannSheet riemann_sheet)
    : k_z(calc_k_z(lm, k_rho, riemann_sheet)),
      Z(dual_solution ? calc_Y(lm, k_z, mode) : calc_Z(lm, k_z, mode)),
      theta(calc_theta(lm, k_z)),
      Gamma_u(calc_Gamma_u(lm, Z, theta)),
      Gamma_d(calc_Gamma_d(lm, Z, theta))
{}

std::vector<cmplx> calc_R(const Internals &d, int n)
{
    std::vector<cmplx> R(4);

    R[0] = d.Gamma_d[n];
    R[1] = d.Gamma_u[n];
    R[2] = d.Gamma_u[n] * d.Gamma_d[n];
    R[3] = R[2];

    return R;
}

std::vector<cmplx> calc_zeta(const LayeredMedium &lm,
                             real z,
                             real z_,
                             int n)
{
    std::vector<cmplx> zeta(4);

    real t = z + z_;

    if (std::isfinite(lm.z[n])) {
        zeta[0] = t - 2 * lm.z[n];
    }

    if (std::isfinite(lm.z[n + 1])) {
        zeta[1] = 2 * lm.z[n + 1] - t;
    }

    if (std::isfinite(lm.d[n])) {
        real t1 = 2 * lm.d[n];
        real t2 = z - z_;
        zeta[2] = t1 + t2;
        zeta[3] = t1 - t2;
    }

    return zeta;
}

cmplx calc_D(const LayeredMedium &lm,
             int n,
             const Internals &d)
{
    using std::complex_literals::operator""i;

    cmplx D = 1.0;
    if (std::isfinite(lm.d[n])) {
        D -= d.Gamma_u[n] * d.Gamma_d[n] * exp(-2.0i * d.theta[n]);
    }

    return D;
}

cmplx V_i_src_wo_prefac(const LayeredMedium &lm,
                        real z,
                        real z_,
                        int n,
                        const Internals &d,
                        bool direct_term)
{
    using std::complex_literals::operator""i;
    auto R = calc_R(d, n);
    auto zeta = calc_zeta(lm, z, z_, n);
    auto D = calc_D(lm, n, d);

    cmplx t1;
    if (direct_term) {
        t1 = std::exp(-1.0i * d.k_z[n] * std::abs(z - z_));
    } else {
        t1 = 0;
    }

    cmplx t2 = 0;
    for (int s = 0; s < 4; s++) {
        t2 += R[s] * std::exp(-1.0i * d.k_z[n] * zeta[s]);
    }
    t2 /= D;

    return (t1 + t2) / 2.0;
}

cmplx V_i_src(const LayeredMedium &lm,
              real z,
              real z_,
              int n,
              const Internals &d,
              bool direct_term)
{
    return d.Z[n] * V_i_src_wo_prefac(lm, z, z_, n, d, direct_term);
}

cmplx I_i_src(const LayeredMedium &lm,
              real z,
              real z_,
              int n,
              const Internals &d,
              bool direct_term)
{
    using std::complex_literals::operator""i;
    auto R = calc_R(d, n);
    auto zeta = calc_zeta(lm, z, z_, n);
    auto D = calc_D(lm, n, d);

    cmplx t1;
    if (direct_term) {
        real dist = z - z_;
        real sign;
        if (dist > 0.0) {
            sign = 1;
        } else if (dist < 0.0) {
            sign = -1;
        } else {
            // Note: set direct term to zero for z == z′.
            // Argument: the function is undefined only for z == z′, i.e., only
            //           at one point => one point contributes nothing to an
           //            integral.
            sign = 0;
        }
        t1 = sign * std::exp(-1.0i * d.k_z[n] * std::abs(dist));
    } else {
        t1 = 0;
    }

    cmplx t2 = 0.0;
    for (int a = 0; a < 4; a++) {
        int s = a + 1;
        t2 += std::pow(-1, s) * R[a] * std::exp(-1.0i * d.k_z[n] * zeta[a]);
    }
    t2 /= D;

    return (t1 - t2) / 2.0;
}

cmplx calc_tau_ud(int n,
                  const std::vector<cmplx> &Gamma_ud,
                  const std::vector<cmplx> &theta)
{
    using std::complex_literals::operator""i;

    return (1.0 + Gamma_ud[n]) * std::exp(-1.0i * theta[n]) /
            (1.0 + Gamma_ud[n] * std::exp(-2.0i * theta[n]));
}

cmplx calc_tau_prod_d(int m, int n, const Internals &d)
{
    cmplx val = 1.0;
    for (int k = m + 1; k <= n - 1; k++) {
        val *= calc_tau_ud(k, d.Gamma_d, d.theta);
    }

    return val;
}

cmplx calc_tau_prod_u(int m, int n, const Internals &d)
{
    cmplx val = 1.0;
    for (int k = n + 1; k <= m - 1; k++) {
        val *= calc_tau_ud(k, d.Gamma_u, d.theta);
    }

    return val;
}

cmplx factor_eq_126(const LayeredMedium &lm,
                    real z,
                    int m,
                    std::function<cmplx(cmplx, cmplx)> pm_operator,
                    const Internals &d)
{
    using std::complex_literals::operator""i;

    cmplx t1 = std::exp(-1.0i * d.k_z[m] * (lm.z[m + 1] - z));

    if (std::isfinite(lm.d[m])) {
        cmplx t2 = 1.0 + d.Gamma_d[m] * std::exp(-2.0i * d.theta[m]);
        cmplx t3 = d.Gamma_d[m] * std::exp(-2.0i * d.k_z[m] * (z - lm.z[m]));
        return t1 / t2 * pm_operator(1, t3);
    } else {
        return t1;
    }
}

cmplx factor_eq_117(const LayeredMedium &lm,
                    real z,
                    int m,
                    std::function<cmplx(cmplx, cmplx)> pm_operator,
                    const Internals &d)
{
    using std::complex_literals::operator""i;
    cmplx t1 = std::exp(-1.0i * d.k_z[m] * (z - lm.z[m]));

    if (std::isfinite(lm.d[m])) {
        cmplx t2 = 1.0 + d.Gamma_u[m] * std::exp(-2.0i * d.theta[m]);
        cmplx t3 = d.Gamma_u[m] * std::exp(-2.0i * d.k_z[m] * (lm.z[m + 1] - z));
        return t1 / t2 * pm_operator(1, t3);
    } else {
        return t1;
    }
}

cmplx T_d(const LayeredMedium &lm,
          real z,
          int m,
          int n,
          std::function<cmplx(cmplx, cmplx)> pm_operator,
          const Internals &d)
{
    return calc_tau_prod_d(m, n, d) *
            factor_eq_126(lm, z, m, pm_operator, d);
}

cmplx T_u(const LayeredMedium &lm,
          real z,
          int m,
          int n,
          std::function<cmplx(cmplx, cmplx)> pm_operator,
          const Internals &d)
{
    return calc_tau_prod_u(m, n, d) *
            factor_eq_117(lm, z, m, pm_operator, d);
}

cmplx V_i_base_wo_prefac(const LayeredMedium &lm,
                         real z,
                         real z_,
                         int m,
                         int n,
                         const Internals &d,
                         bool direct_term)
{
    auto pm_operator = [](cmplx a, cmplx b) { return a + b; };

    cmplx val;
    if (m == n) {
        val = V_i_src_wo_prefac(lm, z, z_, n, d, direct_term);
    } else if (m < n) {
        val = V_i_src_wo_prefac(lm, lm.z[n], z_, n, d, direct_term) *
                T_d(lm, z, m, n, pm_operator, d);
    } else if (m > n) {
        val = V_i_src_wo_prefac(lm, lm.z[n + 1], z_, n, d, direct_term) *
                T_u(lm, z, m, n, pm_operator, d);
    }

    return val;
}

cmplx V_i_base(const LayeredMedium &lm,
               real z,
               real z_,
               int m,
               int n,
               const Internals &d,
               bool direct_term)
{
    return d.Z[n] * V_i_base_wo_prefac(lm, z, z_, m, n, d, direct_term);
}

cmplx I_i_base(const LayeredMedium &lm,
               real z,
               real z_,
               int m,
               int n,
               const Internals &d,
               bool direct_term)
{
    auto pm_operator = [](cmplx a, cmplx b) { return a - b; };

    cmplx val;
    if (m == n) {
        val = I_i_src(lm, z, z_, n, d, direct_term);
    } else if (m < n) {
        val = V_i_src(lm, lm.z[n], z_, n, d, direct_term) *
                T_d(lm, z, m, n, pm_operator, d) / (-d.Z[m]);
    } else if (m > n) {
        val = V_i_src(lm, lm.z[n + 1], z_, n, d, direct_term) *
                T_u(lm, z, m, n, pm_operator, d) / d.Z[m];
    }

    return val;
}

} // namespace utils

} // namespace mthesis::tlgf
