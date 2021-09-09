#include <mthesis.hpp>

#include <gsl/gsl_const_mksa.h>

#include <random>

using namespace mthesis;

VectorR3 rand_point(const FrequencyDomain &fd)
{
    std::random_device dev;
    std::mt19937 gen(dev());
    std::uniform_real_distribution<double> dist_pos(0.0,
                                                    20.0 * fd.lambda_0);

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

    VectorR3 r_ = {0, 0, 0};

    auto r = rand_point(fd);

    std::cout << "r = \n" << r << "\n";
    std::cout << "r_ = \n" << r_ << "\n";

    bool direct_term = false;
    real nu = 0;
    
    auto si = gf::scalar::layered_media::get_sommerfeld_integral(lm, nu, mode,
                                                                 direct_term);

    auto val_ref = si.eval_si_along_sip(r, r_);

    dcim::ThreeLevelV1 dcim_3lv1(si);
    auto val_dcim_3lv1 = dcim_3lv1.get_spatial_gf(r, r_);

    dcim::ThreeLevelV2 dcim_3lv2(si);
    auto val_dcim_3lv2 = dcim_3lv2.get_spatial_gf(r, r_);


    std::cout << "Numerical integration: " << val_ref << "\n";
    std::cout << "DCIM three-level V1:   " << val_dcim_3lv1 << "\n";
    std::cout << "DCIM three-level V2:   " << val_dcim_3lv2 << "\n";
    std::cout << "Error V1:              " <<
                 calc_rel_err_db(val_dcim_3lv1, val_ref) << " dB\n";
    std::cout << "Error V2:              " <<
                 calc_rel_err_db(val_dcim_3lv2, val_ref) << " dB\n";

    return 0;
}
