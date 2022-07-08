#include <mthesis/si/partition_extrapolation.hpp>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/quadrature/gauss.hpp>

#include <cassert>

namespace mthesis::si::pe {

Params::Params() : max_intervals(15), tol(1e-8) {}

void Params::check() const
{
    assert(max_intervals > 0);
    assert(tol > 0.0);
}

real get_first_zero(real nu, real a, real rho)
{
    assert(rho > 0.0);
    assert(a >= 0.0);
    assert(nu >= -1);

    const auto j_table = utils::get_j_table(nu);

    real bessel_arg = a * rho;

    // First search tabulated values.
    real j_found = std::numeric_limits<real>::quiet_NaN(); // Flag.
    for (size_t i = 0; i < j_table.size(); i++) {
        if (j_table[i] > bessel_arg) {
            j_found = j_table[i];
            break;
        }
    }

    // Use McMahon's expansion if not contained in table.
    unsigned m = j_table.size() + 1;
    while (std::isnan(j_found)) {
        auto j = utils::mc_mahon(nu, m);
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

cmplx levin_sidi(std::function<cmplx(real)> f,
                 real rho,
                 real a,
                 Params params)
{
    params.check();
    // Note: Michalski2016a recommends 16-order Gauss-Legendre quadrature
    //	     between consecutive Bessel zeros (see pages 301 and 306).
    constexpr boost::math::quadrature::gauss<real, 15> quad_gl;

    auto xi = utils::get_xi(a, rho, params.max_intervals);
    std::vector<cmplx> A(xi.size() - 1), B(xi.size() - 1);

    cmplx s = 0.0;
    std::vector<cmplx> old(2);

    int k_max = xi.size() - 2;
    for (int k = 0; k <= k_max; k++) {
        cmplx u = quad_gl.integrate(f, xi[k], xi[k + 1]);
        s += u;
        // Detect vanishing function case as levin_sidi_extrap cannot handle
        // the case omega_k = 0 and to avoid failure of the convergence checker.
        if (std::abs(u) < 1e-20 && std::abs(s) < 1e-12) {
            return 0.0;
        }
        cmplx omega = u;
        cmplx val = utils::levin_sidi_extrap(k, s, omega, xi, A, B);
        if (k > 1 && utils::check_converged(val, old, params.tol)) {
            return val;
        }
        old[0] = old[1];
        old[1] = val;
    }

    assert(false);	// Not converged, should never be reached.
}

cmplx mosig_michalski(std::function<cmplx(real)> f,
                      real alpha,
                      real zeta,
                      real rho,
                      real a,
                      Params params)
{
    params.check();
    constexpr boost::math::quadrature::gauss_kronrod<real, 15> quad;

    auto xi = utils::get_xi(a, rho, params.max_intervals);
    std::vector<cmplx> R(xi.size() - 1);

    auto Omega = [&](int k) {
        if (0 == k) {
            return 0.0;
        } else {
            return -std::exp(-(xi[k + 1] - xi[k]) * zeta) *
                    std::pow(xi[k] / xi[k + 1], alpha);
        }
    };

    constexpr real mu = 2.0;
    cmplx s = 0.0;
    std::vector<cmplx> old(2);

    using boost::math::quadrature::gauss_kronrod;
    int k_max = xi.size() - 2;
    for (int k = 0; k <= k_max; k++) {
        cmplx u = quad.integrate(f, xi[k], xi[k + 1]);
        s += u;
        cmplx val = utils::mosig_michalski_extrap(mu, k, s, Omega(k), xi, R);
        if (k > 1 && utils::check_converged(val, old, params.tol)) {
            return val;
        }
        old[0] = old[1];
        old[1] = val;
    }

    assert(false);	// Not converged, should never be reached.
}

namespace utils {

real mc_mahon(real nu, unsigned m)
{
    using namespace std;

    assert(nu >= 0.0);
    real mu = 4.0 * pow(nu, 2.0);
    real a = M_PI * (m + 0.5 * nu - 0.25);

    real j = a -
            (mu - 1) / (8 * a) -
            4 * (mu - 1) * (7 * mu - 31) / (3 * pow(8 * a, 3)) -
            32 * (mu - 1) * (83 * pow(mu, 2) - 982 * mu + 3779) / (15 * pow(8 * a, 5)) -
            64 * (mu - 1) * (6949 * pow(mu, 3) - 153855 * pow(mu, 2) + 1585743 * mu - 6277237) / (105 * pow(8 * a, 7));

    return j;
}

std::vector<real> get_j_table(real nu)
{

    if (0.0 == nu) {
        return std::vector<real>{
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
        return std::vector<real>{
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
        std::vector<real> j_table;
        constexpr unsigned n_roots = 20;
        boost::math::cyl_bessel_j_zero(
                    nu, 1, n_roots, std::back_inserter(j_table));
        return j_table;
    }
}

std::vector<real> get_xi(real a, real rho, unsigned max_intervals)
{
    real q = M_PI / rho;

    std::vector<real> xi(max_intervals + 1);

    xi.front() = a;
    for (size_t i = 1; i < xi.size(); i++) {
        xi[i] = a + i * q;
    }

    return xi;
}

bool check_converged(cmplx val, const std::vector<cmplx> &old, real tol)
{
    std::vector<real> diff(old.size());

    for (size_t i = 0; i < old.size(); i++) {
        diff[i] = std::abs(val - old[i]);
    }

    real max = *std::max_element(diff.begin(), diff.end());

    return max < tol * std::abs(val) ? true : false;
}

cmplx levin_sidi_extrap(int k,
                        cmplx s_k,
                        cmplx omega_k,
                        const std::vector<real> &xi,
                        std::vector<cmplx> &A,
                        std::vector<cmplx> &B)
{
    B[k] = 1.0 / omega_k;
    A[k] = s_k * B[k];
    for (int j = 1; j <= k; j++) {
        real d = 1.0 / xi[k + 1] - 1.0 / xi[k - j + 1];
        A[k - j] = (A[k - j + 1] - A[k - j]) / d;
        B[k - j] = (B[k - j + 1] - B[k - j]) / d;
    }
    return A[0] / B[0];
}

cmplx mosig_michalski_extrap(real mu,
                             int k,
                             cmplx s_k,
                             real Omega_k,
                             const std::vector<real> &xi,
                             std::vector<cmplx> &R)
{
    R[k] = s_k;
    for (int j = 1; j <= k; j++) {
        real d = xi[k - j + 2] - xi[k - j + 1];
        real eta = Omega_k / (1.0 + mu * (j - 1) * d / xi[k - j + 1]);
        R[k - j] = (R[k - j + 1] - eta * R[k - j]) / (1.0 - eta);
    }
    return R[0];
}

} // namespace utils

} // namespace mthesis::si::pe
