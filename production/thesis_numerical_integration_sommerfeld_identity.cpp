#include <mthesis.hpp>

#include <armadillo>
#include <boost/timer/timer.hpp>

#include <filesystem>
#include <cstdio>

using namespace mthesis;

// Note: directory to write file can be passed as single command line argument.
int main(int argc, char** argv)
{
    // Problem specification.
    auto fd = FrequencyDomain(3e9);
    auto vacuum = Vacuum(fd);
    auto lm = FreeSpace(fd);
    EmMode mode = EmMode::TM;	// Does not matter here...
    bool direct_term = true;
    real nu = 0;
    auto si = gf::scalar::layered_media::get_sommerfeld_integral(lm, nu, mode,
                                                                 direct_term);

    VectorR3 r_ = {0, 0, 0};
    r_ *= fd.lambda_0;

    arma::vec x_vals = arma::logspace(-1, 3, 60) * fd.lambda_0;
    arma::vec z_vals = arma::logspace(-1, 3, 60) * fd.lambda_0;

    // Perform computations and measure execution time.
    size_t N = x_vals.size() * z_vals.size();
    std::vector<real> rel_err_db(N);
    std::vector<real> time_s(N);
    boost::timer::auto_cpu_timer timer;
    size_t n = 0;
    for (auto &x : x_vals)
    {
        for (auto &z : z_vals)
        {
            VectorR3 r = {x, 0, z};
            cmplx ref = gf::scalar::free_space::G_0(vacuum, r, r_);

            // Multiple runs for timer.
            int N_runs = 50;
            timer.start();
            for (int i = 0; i < N_runs - 1; i++)
            {
                si.eval_si_along_sip(r, r_);
            }
            // One final time for the result.
            auto num = si.eval_si_along_sip(r, r_);
            timer.stop();

            rel_err_db[n] = calc_rel_err_db(num, ref);
            time_s[n] = std::stod(timer.format(9, "%u"));

            printf("Processing point %4ld of %4ld\n", n + 1, N);
            n++;
        }
    }
    printf("Done.\n");

    printf("Maximum error: %.2f dB\n",
           *std::max_element(rel_err_db.begin(), rel_err_db.end()));

    real time_s_min = *std::min_element(time_s.begin(), time_s.end());

    // Check where to write output.
    FILE *out_target;
    if (2 == argc && std::filesystem::is_directory(argv[1]))
    {
        auto file_fullpath = std::filesystem::path(argv[1]);
        file_fullpath /= "thesis_numerical_integration_sommerfeld_identity.dat";
        out_target = fopen(file_fullpath.c_str(), "w");
    }
    else
    {
        out_target = stdout;
    }

    // Write output.
    fprintf(out_target, "x_by_lambda_0 z_by_lambda_0 rel_err_db rel_time\n");
    n = 0;
    for (auto &x : x_vals)
    {
        for (auto &z : z_vals)
        {
            fprintf(out_target, "%.12f %.12f %.12f %.9f\n",
                    x / fd.lambda_0,
                    z / fd.lambda_0,
                    rel_err_db[n],
                    time_s[n] / time_s_min
                    );
            n++;
        }
        fprintf(out_target, "\n");
    }
    fclose(out_target);

    return 0;
}
