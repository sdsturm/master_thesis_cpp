#ifndef MTHESIS_SOLUTION_DOMAIN_HPP
#define MTHESIS_SOLUTION_DOMAIN_HPP

#include <mthesis/detail/definitions.hpp>

#include <vector>

namespace mthesis
{

    struct FrequencyDomain
    {
        const real f;
        const real omega;
        const real k_0;
        const real lambda_0;

        FrequencyDomain(real f);
    };

    struct Medium
    {
        const FrequencyDomain &fd;
        const cmplx eps_r;
        const cmplx mu_r;
        const cmplx eps;
        const cmplx mu;
        const cmplx k;
        const cmplx eta;

        Medium(const FrequencyDomain &fd, cmplx eps_r, cmplx mu_r);
    };

    cmplx cmplx_permittivity(const FrequencyDomain &fd, real eps_r, real sigma);

    struct Vacuum : public Medium
    {
        Vacuum(const FrequencyDomain &fd);
    };

    enum struct BC
    {
        open,
        PEC,
        PMC
    };

    std::vector<real> layer_thicknesses(const std::vector<real> &in);

    struct LayeredMedium
    {
        const FrequencyDomain &fd;
        const std::vector<real> z;
        const std::vector<real> d;
        const std::vector<Medium> media;
        const BC bottom_bc;
        const BC top_bc;

        LayeredMedium(const FrequencyDomain &fd,
                      const std::vector<real> &z_interfaces,
                      const std::vector<Medium> &media,
                      BC bottom_bc,
                      BC top_bc);

        int identify_layer(real z) const;
    };

    struct HalfSpace : public LayeredMedium
    {
        HalfSpace(const FrequencyDomain &fd, const Medium &ground);
    };

    struct HalfSpacePEC : public LayeredMedium
    {
        HalfSpacePEC(const FrequencyDomain &fd);
    };

    struct HalfSpacePMC : public LayeredMedium
    {
        HalfSpacePMC(const FrequencyDomain &fd);
    };

} // namespace mthesis

#endif
