#include <mthesis/si/partition_extrapolation.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cassert>

namespace mthesis::pe {

void Params::set_max_intervals(int max_intervals)
{
    assert(max_intervals > 0);
    this->max_intervals = max_intervals;
}

void Params::set_tol(double tol)
{
    assert(tol > 0.0);
    this->tol = tol;
}

unsigned Params::get_max_intervals() const
{
    return max_intervals;
}

double Params::get_tol() const
{
    return tol;
}

cmplx levin_sidi(std::function<cmplx(real)> f,
                 real nu,
                 real rho,
                 real a,
                 Params params)
{
    double a_pe_start = utils::get_first_zero(nu, a, rho);

    cmplx gap = utils::integrate_gap(f, a, a_pe_start);

    cmplx tail = utils::levin_sidi_core(f, rho, a_pe_start, params);

    return gap + tail;
}

cmplx mosig_michalski(std::function<cmplx(real)> f,
                      real alpha,
                      real zeta,
                      real nu,
                      real rho,
                      real a,
                      Params params)
{
    double a_pe_start = utils::get_first_zero(nu, a, rho);

    cmplx gap = utils::integrate_gap(f, a, a_pe_start);

    cmplx tail = utils::mosig_michalski_core(f, rho, a, alpha, zeta, params);

    return gap + tail;
}

namespace utils {

double mc_mahon(double nu, unsigned m)
{
    using namespace std;

    assert(nu >= 0.0);
    double mu = 4.0 * pow(nu, 2.0);
    double a = M_PI * (m + 0.5 * nu - 0.25);

    double j = a -
            (mu - 1) / (8 * a) -
            4 * (mu - 1) * (7 * mu - 31) / (3 * pow(8 * a, 3)) -
            32 * (mu - 1) * (83 * pow(mu, 2) - 982 * mu + 3779) / (15 * pow(8 * a, 5)) -
            64 * (mu - 1) * (6949 * pow(mu, 3) - 153855 * pow(mu, 2) + 1585743 * mu - 6277237) / (105 * pow(8 * a, 7));

    return j;
}

std::vector<double> get_j_table(double nu)
{

    if (0.0 == nu) {
        return std::vector<double>{
            2.404825557696,
            5.520078110286,
            8.653727912911,
            11.791534439014,
            14.930917708488,
            18.071063967911,
            21.211636629879,
            24.352471530749,
            27.493479132040,
            30.634606468432,
            33.775820213574,
            36.917098353664,
            40.058425764628,
            43.199791713177,
            46.341188371662,
            49.482609897398,
            52.624051841115,
            55.765510755020,
            58.906983926081,
            62.048469190227};
    } else if (1.0 == nu) {
        return std::vector<double>{
            3.831705970208,
            7.015586669816,
            10.173468135063,
            13.323691936314,
            16.470630050878,
            19.615858510468,
            22.760084380593,
            25.903672087618,
            29.046828534917,
            32.189679910974,
            35.332307550084,
            38.474766234772,
            41.617094212814,
            44.759318997653,
            47.901460887185,
            51.043535183572,
            54.185553641061,
            57.327525437901,
            60.469457845347,
            63.611356698481};
    } else {
        std::vector<double> j_table;
        constexpr unsigned n_roots = 20;
        boost::math::cyl_bessel_j_zero(
                    nu, 1, n_roots, std::back_inserter(j_table));
        return j_table;
    }
}

double get_first_zero(double nu, double a, double rho)
{
    assert(rho > 0.0);
    assert(a >= 0.0);
    assert(nu >= -1);

    const auto j_table = get_j_table(nu);

    double bessel_arg = a * rho;

    // First search tabulated values.
    double j_found = std::numeric_limits<double>::quiet_NaN(); // Flag.
    for (size_t i = 0; i < j_table.size(); i++) {
        if (j_table[i] > bessel_arg) {
            j_found = j_table[i];
            break;
        }
    }

    // Use McMahon's expansion if not contained in table.
    unsigned m = j_table.size();
    while (std::isnan(j_found)) {
        auto j = mc_mahon(nu, m);
        if (j > bessel_arg) {
            j_found = j;
            break;
        } else {
            m++;
        }
    }

    // Normalize zero.
    return j_found / rho;
}

cmplx integrate_gap(std::function<cmplx(real)> f,
                    double a,
                    double a_pe_start)
{
    constexpr boost::math::quadrature::gauss_kronrod<real, 31> quad;
    return quad.integrate(f, a, a_pe_start);
}

std::vector<double> get_xi(double a, double rho, unsigned max_intervals)
{
    double q = M_PI / rho;

    std::vector<double> xi(max_intervals + 1);

    xi.front() = a;
    for (size_t i = 1; i < xi.size(); i++) {
        xi[i] = a + i * q;
    }

    return xi;
}

bool check_converged(cmplx val, const std::vector<cmplx> &old, double tol)
{
    std::vector<double> diff(old.size());

    for (size_t i = 0; i < old.size(); i++) {
        diff[i] = std::abs(val - old[i]);
    }

    double max = *std::max_element(diff.begin(), diff.end());

    return max < tol * std::abs(val) ? true : false;
}

cmplx levin_sidi_extrap(int k,
                        cmplx s_k,
                        cmplx omega_k,
                        const std::vector<double> &xi,
                        std::vector<cmplx> &A,
                        std::vector<cmplx> &B)
{
    B[k] = 1.0 / omega_k;
    A[k] = s_k * B[k];
    for (int j = 1; j <= k; j++) {
        double d = 1.0 / xi[k + 1] - 1.0 / xi[k - j + 1];
        A[k - j] = (A[k - j + 1] - A[k - j]) / d;
        B[k - j] = (B[k - j + 1] - B[k - j]) / d;
    }
    return A[0] / B[0];
}

cmplx levin_sidi_core(std::function<cmplx(real)> f,
                      double rho,
                      double a,
                      Params params)
{
    constexpr boost::math::quadrature::gauss_kronrod<double, 15> quad;

    auto xi = get_xi(a, rho, params.get_max_intervals());
    std::vector<cmplx> A(xi.size() - 1), B(xi.size() - 1);

    cmplx s = 0.0;
    std::vector<cmplx> old(2);

    int k_max = xi.size() - 2;
    for (int k = 0; k <= k_max; k++) {
        cmplx u = quad.integrate(f, xi[k], xi[k + 1]);
        s += u;
        if (std::abs(u) < 1e-100 && std::abs(s) < 1e-20) {
            // Double check if function is identical to zero.
            cmplx u_tot = quad.integrate(f, xi.front(), xi.back());
            assert(std::abs(u_tot) < 1e-16);
            return 0.0;
        }
        cmplx omega = u;
        cmplx val = levin_sidi_extrap(k, s, omega, xi, A, B);
        if (k > 1 && check_converged(val, old, params.get_tol())) {
            return val;
        }
        old[0] = old[1];
        old[1] = val;
    }

    std::cerr << "Levin-Sidi PE quit without reaching convergence.\n";
    return std::numeric_limits<real>::quiet_NaN();
}

cmplx mosig_michalski_extrap(double mu,
                             int k,
                             cmplx s_k,
                             double Omega_k,
                             const std::vector<double> &xi,
                             std::vector<cmplx> &R)
{
    R[k] = s_k;
    for (int j = 1; j <= k; j++) {
        double d = xi[k - j + 2] - xi[k - j + 1];
        double eta = Omega_k / (1.0 + mu * (j - 1) * d / xi[k - j + 1]);
        R[k - j] = (R[k - j + 1] - eta * R[k - j]) / (1.0 - eta);
    }
    return R[0];
}

cmplx mosig_michalski_core(std::function<cmplx(real)> f,
                           double rho,
                           double a,
                           double alpha,
                           double zeta,
                           Params params)
{
    constexpr boost::math::quadrature::gauss_kronrod<double, 15> quad;

    auto xi = get_xi(a, rho, params.get_max_intervals());
    std::vector<cmplx> R(xi.size() - 1);

    auto Omega = [&](int k) {
        if (0 == k) {
            return 0.0;
        } else {
            return -std::exp(-(xi[k + 1] - xi[k]) * zeta) *
                    std::pow(xi[k] / xi[k + 1], alpha);
        }
    };

    constexpr double mu = 2.0;
    cmplx s = 0.0;
    std::vector<cmplx> old(2);

    using boost::math::quadrature::gauss_kronrod;
    int k_max = xi.size() - 2;
    for (int k = 0; k <= k_max; k++) {
        cmplx u = quad.integrate(f, xi[k], xi[k + 1]);
        s += u;
        cmplx val = mosig_michalski_extrap(mu, k, s, Omega(k), xi, R);
        if (k > 1 && check_converged(val, old, params.get_tol())) {
            return val;
        }
        old[0] = old[1];
        old[1] = val;
    }

    std::cerr << "Mosig-Michalski PE quit without reaching convergence.\n";
    return std::numeric_limits<real>::quiet_NaN();
}

} // namespace utils

} // namespace mthesis::pe