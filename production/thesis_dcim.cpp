#include <mthesis.hpp>

#include <armadillo>
#include <boost/timer/timer.hpp>

#include <filesystem>
#include <cstdio>
#include <cassert>

using namespace mthesis;


int main(int argc, char** argv)
{
    // *************************************************************************
    // Specify problem.
    // *************************************************************************

    // Dry ground case from p. 22 in Michalski2016b.
    FrequencyDomain fd(10e6);
    cmplx eps_r = cmplx_permittivity(fd, 3, 0.1e-3);
    Medium ground(fd, eps_r);
    HalfSpace hs(fd, ground);

    bool direct_term = false;	// Direct term can be integrated analytically.
    real nu = 0;

    // Source location.
    real z_ = 1.0 * fd.lambda_0;
    VectorR3 r_ = {0, 0, z_};

    // Observation locations.
    arma::vec z_vals = arma::linspace(0, 25, 70) * fd.lambda_0;
    arma::vec rho_vals= arma::logspace(-1, 3, 70) * fd.lambda_0;

    // Setup of Sommerfeld integrals and routines.
    using mthesis::si::axial_transmission::eval_si_along_sip;
    using mthesis::gf::scalar::layered_media::get_sommerfeld_integral;

    auto si_tm = get_sommerfeld_integral(hs, nu, EmMode::TM, direct_term);
    auto si_te = get_sommerfeld_integral(hs, nu, EmMode::TE, direct_term);

    auto dcim_tm = dcim::ThreeLevelV2(si_tm);
    auto dcim_te = dcim::ThreeLevelV2(si_te);



    // *************************************************************************
    // Set up files.
    // *************************************************************************

    // Check output directory.
    assert(2 == argc);
    assert(std::filesystem::is_directory(argv[1]));
    auto target_dir = std::filesystem::path(argv[1]);

    // Set up file for error results.
    auto err_file= target_dir;
    err_file/= "thesis_dcim_error.dat";
    FILE *err_file_ptr = fopen(err_file.c_str(), "w");

    // Set up file for number of images.
    auto imgs_file = target_dir;
    imgs_file /= "thesis_dcim_images.dat";
    FILE *imgs_file_ptr = fopen(imgs_file.c_str(), "w");

    fprintf(err_file_ptr, "z_by_lambda_0 rho_by_lambda_0 rel_err_db_tm rel_err_db_te abs_ref_tm abs_ref_te\n");
    fprintf(imgs_file_ptr, "z_by_lambda_0 n_images_tm n_images_te\n");


    // *************************************************************************
    // Computations and in-place writing.
    // *************************************************************************
    int ctr = 0;
    int N_pts = z_vals.n_elem * rho_vals.n_elem;
    boost::timer::auto_cpu_timer dcim_timer;
    boost::timer::auto_cpu_timer ref_timer;
    bool ref_timer_on = false;
    real max_err_tm = -std::numeric_limits<real>::infinity();
    real max_err_te = -std::numeric_limits<real>::infinity();
    for (const auto &z : z_vals) {

        dcim_timer.resume();
        // Compute complex images (independent of rho)
        auto ce_vecs_levels_tm = dcim_tm.get_exponentials(z, z_);
        auto ce_vecs_levels_te = dcim_te.get_exponentials(z, z_);
        dcim_timer.stop();

        // Count images and print to respective file.
        int N_imgs_tm = 0;
        int N_imgs_te = 0;
        for (const auto &ce_vec : ce_vecs_levels_tm) {
            N_imgs_tm += ce_vec.size();
        }
        for (const auto &ce_vec : ce_vecs_levels_te) {
            N_imgs_te += ce_vec.size();
        }
        fprintf(imgs_file_ptr, "%.6f %d %d\n",
                z / fd.lambda_0,
                N_imgs_tm,
                N_imgs_te);

        // Evaluate for all rho.
        for (const auto &rho : rho_vals) {

            if (!ref_timer_on) { // To avoid measuring first DCIM run.
                ref_timer_on = true;
                ref_timer.start();
            }
            ref_timer.resume();
            auto val_ref_tm = eval_si_along_sip(si_tm, rho, z, z_);
            auto val_ref_te = eval_si_along_sip(si_te, rho, z, z_);
            ref_timer.stop();

            dcim_timer.resume();
            auto val_dcim_tm = dcim_tm.get_spatial_gf(ce_vecs_levels_tm, rho);
            auto val_dcim_te = dcim_te.get_spatial_gf(ce_vecs_levels_te, rho);
            dcim_timer.stop();

            real rel_err_db_tm = calc_rel_err_db(val_dcim_tm, val_ref_tm);
            real rel_err_db_te = calc_rel_err_db(val_dcim_te, val_ref_te);

            // Set to NaN in case of failure to ignore this points in the plot.
            if (rel_err_db_tm > max_err_tm) {
                max_err_tm = rel_err_db_tm;
            }
            if (rel_err_db_te > max_err_te) {
                max_err_te = rel_err_db_te;
            }

            fprintf(err_file_ptr, "%.6e %.6e %.6e %.6e %.6e %.6e\n",
                    z / fd.lambda_0,
                    rho / fd.lambda_0,
                    rel_err_db_tm,
                    rel_err_db_te,
                    abs(val_ref_tm),
                    abs(val_ref_te));

            printf("Processing point %4d of %4d\n", ctr++, N_pts);
        }
        fprintf(err_file_ptr, "\n");

    } // z-loop


    // Clean up.
    fclose(err_file_ptr);
    fclose(imgs_file_ptr);

    // Print execution time ratio and maximum errors.
    real time_ratio = std::stod(dcim_timer.format(9, "%u")) /
            std::stod(ref_timer.format(9, "%u"));
    printf("\n");
    printf("Execution times:\n");
    printf("t_dcim / t_ref = %f\n", time_ratio);
    printf("t_ref / t_dcim = %f\n", 1.0 / time_ratio);
    printf("\n");
    printf("Maximum errors:\n");
    printf("TM: %.2f dB\n", max_err_tm);
    printf("TE: %.2f dB\n", max_err_te);
    printf("\n");

    return 0;
}
