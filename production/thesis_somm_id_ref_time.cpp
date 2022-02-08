#include <mthesis.hpp>

#include <armadillo>

#include <chrono>
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
    size_t n = 0;
    using mthesis::si::axial_transmission::eval_si_along_sip;
    using t_unit = std::chrono::nanoseconds;
    constexpr real t_scale = 1e9;
    for (auto &x : x_vals)
    {
        for (auto &z : z_vals)
        {
            VectorR3 r = {x, 0, z};
            cmplx ref = gf::scalar::free_space::G_0(vacuum, r, r_);

            // Multiple runs for timer.
            int N_runs = 50;
            const auto t1 = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < N_runs - 1; i++)
            {
                eval_si_along_sip(si, r, r_);
            }
            auto num = eval_si_along_sip(si, r, r_);
            const auto t2 = std::chrono::high_resolution_clock::now();
            const auto t = std::chrono::duration_cast<t_unit>(t2 - t1);

            rel_err_db[n] = calc_rel_err_db(num, ref);
            time_s[n] = t.count() / t_scale;

            printf("Processing point %4ld of %4ld, err = %.2f db\n",
                   n + 1,
                   N,
                   rel_err_db[n]);

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
        file_fullpath /= "thesis_somm_id_ref_time.dat";
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
            fprintf(out_target, "%.12e %.12e %.12e %.9e\n",
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
