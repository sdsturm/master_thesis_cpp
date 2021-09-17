#include <mthesis/tlgf.hpp>

namespace mthesis::tlgf {

TLGFParams::TLGFParams(const LayeredMedium &lm,
                       EmMode mode,
                       bool direct_term)
    : lm(lm), mode(mode), direct_term(direct_term)
{}

cmplx V_i(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet)
{
    bool dual_solution = false;
    utils::LayerCoords c(p.lm, z, z_);
    utils::Internals d(p, k_rho, dual_solution, sheet);
    return utils::V_i_base(p, c, d);
}

cmplx V_i(cmplx k_rho, real z, real z_, TLGFParams p)
{
    RiemannSheet proper_sheet = RiemannSheet::I;
    return V_i(k_rho, z, z_, p, proper_sheet);
}

cmplx I_i(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet)
{
    bool dual_solution = false;
    utils::LayerCoords c(p.lm, z, z_);
    utils::Internals d(p, k_rho, dual_solution, sheet);
    return utils::I_i_base(p, c, d);
}

cmplx I_i(cmplx k_rho, real z, real z_, TLGFParams p)
{
    RiemannSheet proper_sheet = RiemannSheet::I;
    return I_i(k_rho, z, z_, p, proper_sheet);
}

cmplx I_v(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet)
{
    bool dual_solution = true;
    utils::LayerCoords c(p.lm, z, z_);
    utils::Internals d(p, k_rho, dual_solution, sheet);
    return utils::V_i_base(p, c, d);
}

cmplx I_v(cmplx k_rho, real z, real z_, TLGFParams p)
{
    RiemannSheet proper_sheet = RiemannSheet::I;
    return I_v(k_rho, z, z_, p, proper_sheet);
}

cmplx V_v(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet)
{
    bool dual_solution = true;
    utils::LayerCoords c(p.lm, z, z_);
    utils::Internals d(p, k_rho, dual_solution, sheet);
    return utils::I_i_base(p, c, d);
}

cmplx V_v(cmplx k_rho, real z, real z_, TLGFParams p)
{
    RiemannSheet proper_sheet = RiemannSheet::I;
    return V_v(k_rho, z, z_, p, proper_sheet);
}

cmplx generic_sgf(cmplx k_rho, real z, real z_, TLGFParams p, RiemannSheet sheet)
{
    using std::complex_literals::operator""i;

    bool dual_solution = false;
    utils::LayerCoords c(p.lm, z, z_);
    utils::Internals d(p, k_rho, dual_solution, sheet);

    // Note: this function is not reciprocal if z and z_ lie in layers with
    // 	     different material paramters.
    //       For the case m == n (z and z_ in the same layer) this function
    //       fulfills reciprocity.
    return utils::V_i_base_wo_prefac(p, c, d) / (1.0i * d.k_z[c.n]);
}

cmplx generic_sgf(cmplx k_rho, real z, real z_, TLGFParams p)
{
    RiemannSheet proper_sheet = RiemannSheet::I;
    return generic_sgf(k_rho, z, z_, p, proper_sheet);
}

// *****************************************************************************

namespace utils {

void calc_k_z_select_sheet(std::vector<cmplx> &k_z, RiemannSheet sheet)
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

std::vector<cmplx> calc_k_z(const TLGFParams &p,
                            cmplx k_rho,
                            RiemannSheet sheet)
{
    std::vector<cmplx> k_z(p.lm.media.size());
    cmplx k_rho_square = std::pow(k_rho, 2.0);

    for (size_t n = 0; n < p.lm.media.size(); n++) {
        k_z[n] = std::sqrt(std::pow(p.lm.media[n].k, 2.0) - k_rho_square);
    }

    calc_k_z_select_sheet(k_z, sheet);

    return k_z;
}

std::vector<cmplx> calc_Z(const TLGFParams &p, const std::vector<cmplx> &k_z)
{
    std::vector<cmplx> Z(p.lm.media.size());

    switch (p.mode) {
    case EmMode::TM:
        for (size_t n = 0; n < p.lm.media.size(); n++) {
            Z[n] = k_z[n] / (p.lm.fd.omega * p.lm.media[n].eps);
        }
        break;
    case EmMode::TE:
        for (size_t n = 0; n < p.lm.media.size(); n++) {
            Z[n] = (p.lm.fd.omega * p.lm.media[n].mu) / k_z[n];
        }
        break;
    }

    return Z;
}

std::vector<cmplx> calc_Y(const TLGFParams &p, const std::vector<cmplx> &k_z)
{
    auto Y = calc_Z(p, k_z);
    for (auto &elem : Y) {
        elem = 1.0 / elem;
    }

    return Y;
}

std::vector<cmplx> calc_theta(const TLGFParams &p,
                              const std::vector<cmplx> &k_z)
{
    std::vector<cmplx> theta(p.lm.media.size());

    for (size_t n = 0; n < p.lm.media.size(); n++) {
        theta[n] = k_z[n] * p.lm.d[n];
    }

    return theta;
}

cmplx calc_Gamma_fresnel(const std::vector<cmplx> &Z, int i, int j)
{
    return (Z[i] - Z[j]) / (Z[i] + Z[j]);
}

std::vector<cmplx> calc_Gamma_u(const TLGFParams &p,
                                const std::vector<cmplx> &Z,
                                const std::vector<cmplx> &theta)
{
    std::vector<cmplx> Gamma_u(p.lm.media.size());

    switch (p.lm.top_bc) {
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
    for (int n = p.lm.media.size() - 2; n >= 0; n--) {
        auto t1 = calc_Gamma_fresnel(Z, n + 1, n);
        if (std::isfinite(p.lm.d[n + 1])) {
            cmplx t2 = Gamma_u[n + 1] * std::exp(-2.0i * theta[n + 1]);
            Gamma_u[n] = (t1 + t2) / (1.0 + t1 * t2);
        } else {
            Gamma_u[n] = t1;
        }
    }

    return Gamma_u;
}

std::vector<cmplx> calc_Gamma_d(const TLGFParams &p,
                                const std::vector<cmplx> &Z,
                                const std::vector<cmplx> &theta)
{
    std::vector<cmplx> Gamma_d(p.lm.media.size());

    switch (p.lm.bottom_bc) {
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
    for (size_t n = 1; n < p.lm.media.size(); n++) {
        auto t1 = calc_Gamma_fresnel(Z, n - 1, n);
        if (std::isfinite(p.lm.d[n - 1])) {
            cmplx t2 = Gamma_d[n - 1] * std::exp(-2.0i * theta[n - 1]);
            Gamma_d[n] = (t1 + t2) / (1.0 + t1 * t2);
        } else {
            Gamma_d[n] = t1;
        }
    }

    return Gamma_d;
}

Internals::Internals(const TLGFParams &p,
                     cmplx k_rho,
                     bool dual_solution,
                     RiemannSheet sheet)
    : k_z(calc_k_z(p, k_rho, sheet)),
      Z(dual_solution ? calc_Y(p, k_z) : calc_Z(p, k_z)),
      theta(calc_theta(p, k_z)),
      Gamma_u(calc_Gamma_u(p , Z, theta)),
      Gamma_d(calc_Gamma_d(p, Z, theta))
{}

LayerCoords::LayerCoords(const LayeredMedium &lm, real z, real z_)
    : z(z), z_(z_), m(lm.identify_layer(z)), n(lm.identify_layer(z_))
{}

std::vector<cmplx> calc_R(const Internals &d, const LayerCoords &c)
{
    std::vector<cmplx> R(4);

    R[0] = d.Gamma_d[c.n];
    R[1] = d.Gamma_u[c.n];
    R[2] = d.Gamma_u[c.n] * d.Gamma_d[c.n];
    R[3] = R[2];

    return R;
}

std::vector<cmplx> calc_zeta(const TLGFParams &p, const LayerCoords &c, real z)
{
    std::vector<cmplx> zeta(4);

    real t = z + c.z_;

    if (std::isfinite(p.lm.z[c.n])) {
        zeta[0] = t - 2 * p.lm.z[c.n];
    }

    if (std::isfinite(p.lm.z[c.n + 1])) {
        zeta[1] = 2 * p.lm.z[c.n + 1] - t;
    }

    if (std::isfinite(p.lm.d[c.n])) {
        real t1 = 2 * p.lm.d[c.n];
        real t2 = z - c.z_;
        zeta[2] = t1 + t2;
        zeta[3] = t1 - t2;
    }

    return zeta;
}

cmplx calc_D(const TLGFParams &p, const LayerCoords &c, const Internals &d)
{
    using std::complex_literals::operator""i;

    cmplx D = 1.0;
    if (std::isfinite(p.lm.d[c.n])) {
        D -= d.Gamma_u[c.n] * d.Gamma_d[c.n] * exp(-2.0i * d.theta[c.n]);
    }

    return D;
}

cmplx V_i_src_wo_prefac(const TLGFParams &p,
                        const LayerCoords &c,
                        real z,
                        const Internals &d)
{
    using std::complex_literals::operator""i;
    auto R = calc_R(d, c);
    auto zeta = calc_zeta(p, c, z);
    auto D = calc_D(p, c, d);

    cmplx t1;
    if (p.direct_term) {
        t1 = std::exp(-1.0i * d.k_z[c.n] * std::abs(z - c.z_));
    } else {
        t1 = 0;
    }

    cmplx t2 = 0;
    for (int s = 0; s < 4; s++) {
        t2 += R[s] * std::exp(-1.0i * d.k_z[c.n] * zeta[s]);
    }
    t2 /= D;

    return (t1 + t2) / 2.0;
}

cmplx V_i_src(const TLGFParams &p,
              const LayerCoords &c,
              real z,
              const Internals &d)
{
    return d.Z[c.n] * V_i_src_wo_prefac(p, c, z, d);
}

cmplx I_i_src(const TLGFParams &p,
              const LayerCoords &c,
              real z,
              const Internals &d)
{
    using std::complex_literals::operator""i;
    auto R = calc_R(d, c);
    auto zeta = calc_zeta(p, c, z);
    auto D = calc_D(p, c, d);

    cmplx t1;
    if (p.direct_term) {
        real dist = z - c.z_;
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
        t1 = sign * std::exp(-1.0i * d.k_z[c.n] * std::abs(dist));
    } else {
        t1 = 0;
    }

    cmplx t2 = 0.0;
    for (int a = 0; a < 4; a++) {
        int s = a + 1;
        t2 += std::pow(-1, s) * R[a] * std::exp(-1.0i * d.k_z[c.n] * zeta[a]);
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

cmplx calc_tau_prod_d(const LayerCoords &c, const Internals &d)
{
    cmplx val = 1.0;
    for (int k = c.m + 1; k <= c.n - 1; k++) {
        val *= calc_tau_ud(k, d.Gamma_d, d.theta);
    }

    return val;
}

cmplx calc_tau_prod_u(const LayerCoords &c, const Internals &d)
{
    cmplx val = 1.0;
    for (int k = c.n + 1; k <= c.m - 1; k++) {
        val *= calc_tau_ud(k, d.Gamma_u, d.theta);
    }

    return val;
}

cmplx factor_eq_126(const TLGFParams &p,
                    const LayerCoords &c,
                    std::function<cmplx(cmplx, cmplx)> pm_operator,
                    const Internals &d)
{
    using std::complex_literals::operator""i;

    cmplx t1 = std::exp(-1.0i * d.k_z[c.m] * (p.lm.z[c.m + 1] - c.z));

    if (std::isfinite(p.lm.d[c.m])) {
        cmplx t2 = 1.0 + d.Gamma_d[c.m] * std::exp(-2.0i * d.theta[c.m]);
        cmplx t3 = d.Gamma_d[c.m] *
                std::exp(-2.0i * d.k_z[c.m] * (c.z - p.lm.z[c.m]));
        return t1 / t2 * pm_operator(1, t3);
    } else {
        return t1;
    }
}

cmplx factor_eq_117(const TLGFParams &p,
                    const LayerCoords &c,
                    std::function<cmplx(cmplx, cmplx)> pm_operator,
                    const Internals &d)
{
    using std::complex_literals::operator""i;
    cmplx t1 = std::exp(-1.0i * d.k_z[c.m] * (c.z - p.lm.z[c.m]));

    if (std::isfinite(p.lm.d[c.m])) {
        cmplx t2 = 1.0 + d.Gamma_u[c.m] * std::exp(-2.0i * d.theta[c.m]);
        cmplx t3 = d.Gamma_u[c.m] *
                std::exp(-2.0i * d.k_z[c.m] * (p.lm.z[c.m + 1] - c.z));
        return t1 / t2 * pm_operator(1, t3);
    } else {
        return t1;
    }
}

cmplx T_d(const TLGFParams &p,
          const LayerCoords &c,
          std::function<cmplx(cmplx, cmplx)> pm_operator,
          const Internals &d)
{
    return calc_tau_prod_d(c, d) * factor_eq_126(p, c, pm_operator, d);
}

cmplx T_u(const TLGFParams &p,
          const LayerCoords &c,
          std::function<cmplx(cmplx, cmplx)> pm_operator,
          const Internals &d)
{
    return calc_tau_prod_u(c, d) * factor_eq_117(p, c, pm_operator, d);
}

cmplx V_i_base_wo_prefac(const TLGFParams &p,
                         const LayerCoords &c,
                         const Internals &d)
{
    auto pm_operator = [](cmplx a, cmplx b) { return a + b; };

    cmplx val;
    if (c.m == c.n) {
        val = V_i_src_wo_prefac(p, c, c.z, d);
    } else if (c.m < c.n) {
        val = V_i_src_wo_prefac(p, c, p.lm.z[c.n], d) *
                T_d(p, c, pm_operator, d);
    } else if (c.m > c.n) {
        val = V_i_src_wo_prefac(p, c, p.lm.z[c.n + 1], d) *
                T_u(p, c, pm_operator, d);
    }

    return val;
}

cmplx V_i_base(const TLGFParams &p,
               const LayerCoords &c,
               const Internals &d)
{
    return d.Z[c.n] * V_i_base_wo_prefac(p, c, d);
}

cmplx I_i_base(const TLGFParams &p,
               const LayerCoords &c,
               const Internals &d)
{
    auto pm_operator = [](cmplx a, cmplx b) { return a - b; };

    cmplx val;
    if (c.m == c.n) {
        val = I_i_src(p, c, c.z, d);
    } else if (c.m < c.n) {
        val = V_i_src(p, c, p.lm.z[c.n], d) *
                T_d(p, c, pm_operator, d) / (-d.Z[c.m]);
    } else if (c.m > c.n) {
        val = V_i_src(p, c, p.lm.z[c.n + 1], d) *
                T_u(p, c, pm_operator, d) / d.Z[c.m];
    }

    return val;
}

}

// namespace utils

} // namespace mthesis::tlgf
