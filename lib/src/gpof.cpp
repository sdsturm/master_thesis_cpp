#include <mthesis/detail/gpof.hpp>

#include <armadillo>

#include <cassert>

namespace mthesis
{
    std::vector<ComplexExponential> gpof(const std::vector<cmplx> &y,
                                         real d_t,
                                         const GPOFParams params)
    {
        assert(d_t > 0.0);
        int N = y.size();

        // Set pencil parameter.
        int L;
        if (params.L > 0)
        {
            L = params.L;
        }
        else
        {
            L = N / 2 - 1;
        }

        // Build matrix Y.
        arma::cx_mat Y(N - L, L + 1);
        for (arma::uword i = 0; i < Y.n_rows; i++)
            for (arma::uword j = 0; j < Y.n_cols; j++)
                Y(i, j) = y[i + j];

        // Perform SVD.
        arma::cx_mat U;
        arma::vec s;
        arma::cx_mat V;
        arma::svd_econ(U, s, V, Y, "right");

        // Determine model order.
        int M;
        if (params.M < 0)
        {
            M = -1;
            for (arma::uword i = 1; i < s.n_elem; i++)
            {
                if (std::abs(s(i) / s(1)) < params.tol)
                {
                    M = i;
                    break;
                }
            }
            if (-1 == M)
            {
                M = std::min(params.M_max, L - 1);
            }
        }
        else
        {
            M = params.M;
        }
        assert(M <= L);
        assert(L <= N - M);

        // Compute eigenvalues.
        arma::cx_mat Y_1 =
            V(arma::span(1, V.n_rows - 1), arma::span(0, M - 1)).t() *
            V(arma::span(0, V.n_rows - 2), arma::span(0, M - 1));

        arma::cx_mat Y_2 =
            V(arma::span(0, V.n_rows - 2), arma::span(0, M - 1)).t() *
            V(arma::span(0, V.n_rows - 2), arma::span(0, M - 1));

        arma::cx_vec z = arma::eig_pair(Y_1, Y_2);

        // Compute exponents (see Appendix I in Sarkar1995 not (21) and (22)).
        arma::cx_mat Z(N, z.n_elem);
        for (int r = 0; r < N; r++)
            for (int c = 0; c < static_cast<int>(z.n_elem); c++)
                Z(r, c) = std::pow(z(c), r);

        arma::cx_vec rhs(y);
        arma::cx_vec b = arma::solve(Z, rhs);

        // Postprocessing.
        assert(b.n_elem == z.n_elem);

        std::vector<ComplexExponential> ans;
        for (arma::uword i = 0; i < z.n_elem; i++)
            ans.emplace_back(b(i), std::log(z(i)) / d_t);

        return ans;
    }

    std::vector<cmplx> reconstruct_signal(std::vector<ComplexExponential> &ce,
                                          double d_t,
                                          unsigned N)
    {
        std::vector<cmplx> y(N);

        for (size_t n = 0; n < N; n++)
        {
            for (const auto &elem : ce)
            {
                double t = n * d_t;
                y[n] += elem.amplitude * std::exp(elem.exponent * t);
            }
        }
        return y;
    }

} // namespace mthesis
