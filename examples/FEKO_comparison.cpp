#include <mthesis.hpp>

#include <armadillo>

namespace mthesis {

void get_feko_data(std::string filename,
                   std::vector<VectorR3> &obs_pts,
                   std::vector<VectorC3> &field_pts)
{
    std::ifstream feko_file;
    feko_file.open(filename);
    if (!feko_file.is_open()) {
        throw  std::runtime_error("Could not open file.");
    }

    // Get rid of header
    std::string line;
    do {
        std::getline(feko_file, line);
    } while (line.find("#No. of Header Lines: ") == std::string::npos);
    std::getline(feko_file, line);

    // Parse all lines
    std::vector<double> vals(9);
    while (std::getline(feko_file, line)) {
        std::stringstream ss(line);
        int col_index = 0;
        while (!ss.eof()) {
            while(ss.peek() == ' ') {
                ss.ignore();
            }
            ss >> vals[col_index];
            col_index++;
        }
        VectorR3 r = {vals[0], vals[1], vals[2]};
        VectorC3 F = {cmplx(vals[3], vals[4]), cmplx(vals[5], vals[6]), cmplx(vals[7], vals[8])};
        obs_pts.push_back(r);
        field_pts.push_back(F);
    }

    feko_file.close();
}

} // namespace mthesis

using namespace mthesis;

int main()
{
    // Import FEKO data
    std::vector<VectorR3> r_vals;
    std::vector<VectorC3> E_vals_FEKO;
//    get_feko_data("/home/sistu/HFT/research/codes/cpp-libs/HSGF/examples/dry_ground_VED_h_1_lambda0_cartesian_sampling.efe", r_vals, E_vals_FEKO);
    get_feko_data("/home/sistu/HFT/research/codes/cpp-libs/HSGF/examples/dry_ground_VED_h_1_lambda0_cartesian_sampling.hfe", r_vals, E_vals_FEKO);

    // Problem specification
    FrequencyDomain fd(10e6);
    auto eps_r = cmplx_permittivity(fd, 3, 0.1e-3);
    Medium ground(fd, eps_r);
    HalfSpace hs(fd, ground);
    VectorR3 r_ = {0.0, 0.0, 1.0 * fd.lambda_0};
    VectorR3 dipole_dir = {0.0, 0.0, 1.0};
    cmplx dipole_Il = 1.0;

    // Compute own solution
    std::vector<VectorC3> E_vals_sistu(r_vals.size());
//#pragma omp parallel for
    for (size_t n = 0; n < r_vals.size(); ++n) {
//        DyadC3 G_EJ = gf::dyadic::layered_media::G_EJ(hs, r_vals[n], r_);
        DyadC3 G_EJ = gf::dyadic::layered_media::G_HJ(hs, r_vals[n], r_);
        E_vals_sistu[n](0) = G_EJ(0, 2);
        E_vals_sistu[n](1) = G_EJ(1, 2);
        E_vals_sistu[n](2) = G_EJ(2, 2);
    }

    // Print results
    real max_err = 0;
    for (size_t n = 0; n < r_vals.size(); ++n) {
        std::cout << "r =\n" << r_vals[n];
        std::cout << "FEKO:\n" << E_vals_FEKO[n];
        std::cout << "SiStu:\n" << E_vals_sistu[n];
        real err = arma::norm(E_vals_sistu[n] - E_vals_FEKO[n]) / arma::norm(E_vals_FEKO[n]);
        if (err > max_err) max_err = err;
        std::cout << "max abs diff: " << 20 * log10(err) << "dB\n";
        std::cout << "================================================================================\n";
    }

    std::cout << "\n" << "Max relative error: " << 20 * log10(max_err) << "dB\n";

    return 0;
}
