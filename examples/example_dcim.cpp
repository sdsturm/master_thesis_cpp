#include <mthesis.hpp>

#include <gsl/gsl_const_mksa.h>

#include <random>

using namespace mthesis;

VectorR3 rand_point(real dist)
{
    std::random_device dev;
    std::mt19937 gen(dev());
    std::uniform_real_distribution<double> dist_pos(0.0, dist);

    VectorR3 r = {dist_pos(gen), dist_pos(gen), dist_pos(gen)};

    return r;
}

int main()
{
    using std::complex_literals::operator""i;

    auto fd = FrequencyDomain(GSL_CONST_MKSA_SPEED_OF_LIGHT / 1.0);
    auto ground = Medium(fd, 3.0, 1.0);
    auto lm = HalfSpace(fd, ground);
    auto mode = EmMode::TM;

     auto r = rand_point(5 * fd.lambda_0);
     VectorR3 r_ = {0, 0, 0};
//     auto r_ = rand_point(5 * fd.lambda_0);

    std::cout << "r = \n" << r << "\n";
    std::cout << "r_ = \n" << r_ << "\n";

    bool direct_term = false;
    real nu = 0;
    
    auto si = gf::scalar::layered_media::get_sommerfeld_integral(lm, nu, mode,
                                                                 direct_term);

    auto lmc = LayeredMediumCoords(r, r_);
    auto val_ref = si.eval_si_along_sip(r, r_);

    dcim::ThreeLevelV2 my_dcim(si);
    auto ce_vecs_levels = my_dcim.get_exponentials(lmc.z, lmc.z_);
    auto val_dcim = my_dcim.get_spatial_gf(ce_vecs_levels, lmc.rho);

    std::cout << "reference: " << val_ref << "\n";
    std::cout << "dcim:      " << val_dcim << "\n";
    std::cout << "Error:     " << calc_rel_err_db(val_dcim, val_ref) << " dB\n";

    return 0;
}
