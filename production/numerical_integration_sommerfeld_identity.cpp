#include <mthesis.hpp>

#include <armadillo>
#include <boost/timer/timer.hpp>

#include <filesystem>
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

    auto data_file_path = std::filesystem::path(MTHESIS_RESULTS_DIR);
    data_file_path /= "numerical_integration_sommerfeld_identity.dat";
    FILE *data_file;
    data_file = fopen(data_file_path.c_str(), "w");
    if (!data_file)
    {
        std::cerr << "Unable to open TEX data output file.\n";
        return -1;
    }
    fprintf(data_file, "x_by_lambda_0 z_by_lambda_0 rel_err_db time\n");
    size_t n = 1;
    size_t N = x_vals.size() * z_vals.size();
    boost::timer::auto_cpu_timer timer;
    for (auto &x : x_vals)
    {
        for (auto &z : z_vals)
        {
            VectorR3 r = {x, 0, z};
            auto ref = sgf::free_space(vacuum, r, r_);

            // Multiple runs for timer.
            int N_runs = 50;
            timer.start();
            for (int i = 0; i < N_runs - 1; i++)
                sgf::lm_generic_spatial(lm, r, r_, mode, direct_term);
            // One final time for the result.
            auto num = sgf::lm_generic_spatial(lm, r, r_, mode, direct_term);
            timer.stop();

            auto rel_err_db = calc_rel_err_db(num, ref);
            auto t_str = timer.format(9, "%u");

            fprintf(data_file, "%.12f %.12f %.12f %s\n",
                    x / fd.lambda_0,
                    z / fd.lambda_0,
                    rel_err_db,
                    t_str.c_str());

            printf("Processing point %4ld of %4ld\n", n, N);
            n++;
        }
        fprintf(data_file, "\n");
    }
    printf("Done.\n");
    fclose(data_file);

    return 0;
}
