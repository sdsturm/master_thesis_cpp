#include <mthesis/mthesis.hpp>

#include <armadillo>

#include <cstdio>

using namespace mthesis;

int main()
{
    auto fd = FrequencyDomain(3e9);
    auto vacuum = Vacuum(fd);
    auto lm = FreeSpace(fd);
    EmMode mode = EmMode::TM;
    bool direct_term = true;

    VectorR3 r_ = {0, 0, 0};
    r_ *= fd.lambda_0;

    arma::vec x_vals = arma::logspace(-1, 3, 60) * fd.lambda_0;
    arma::vec z_vals = arma::logspace(-1, 3, 60) * fd.lambda_0;

    real max_rel_err_db = -INFINITY;

    FILE *pFile;
    pFile = fopen("/home/sistu/Documents/test/data.txt", "w");
    fprintf(pFile, "x_by_lambda_0  z_by_lambda_0  rel_err_db \n");

    for (auto &x : x_vals)
    {
        for (auto &z : z_vals)
        {
            VectorR3 r = {x, 0, z};
            auto num = sgf::lm_generic_spatial(lm, r, r_, mode, direct_term);
            auto ref = sgf::free_space(vacuum, r, r_);

            auto rel_err_db = calc_rel_err_db(num, ref);
            printf("x = %.2e   z = %.2e   err = %.2f dB\n", x, z, rel_err_db);

            fprintf(pFile, "%.12f  %.12f  %.12f  \n",
                    x / fd.lambda_0,
                    z / fd.lambda_0,
                    rel_err_db);
            
            if (rel_err_db > max_rel_err_db)
                max_rel_err_db = rel_err_db;
        }
        fprintf(pFile, "\n");
    }
    fclose(pFile);

    printf("\n\nMaximum error: %.2f dB\n", max_rel_err_db);

    return 0;
}
