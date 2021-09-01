#include <mthesis/detail/solution_domain.hpp>

#include <gsl/gsl_const_mksa.h>

#include <cassert>

namespace mthesis
{
    FrequencyDomain::FrequencyDomain(real f)
        : f(f),
          omega(2.0 * M_PI * f),
          k_0(omega / GSL_CONST_MKSA_SPEED_OF_LIGHT),
          lambda_0(2.0 * M_PI / k_0)
    {
    }

    Medium::Medium(const FrequencyDomain &fd, cmplx eps_r, cmplx mu_r)
        : fd(fd),
          eps_r(eps_r),
          mu_r(mu_r),
          eps(GSL_CONST_MKSA_VACUUM_PERMITTIVITY * eps_r),
          mu(GSL_CONST_MKSA_VACUUM_PERMEABILITY * mu_r),
          k(fd.omega * std::sqrt(eps * mu)),
          eta(sqrt(mu / eps))
    {
    }

    cmplx cmplx_permittivity(const FrequencyDomain &fd, real eps_r, real sigma)
    {
        return cmplx(eps_r,
                     -sigma / (fd.omega * GSL_CONST_MKSA_VACUUM_PERMITTIVITY));
    }

    Vacuum::Vacuum(const FrequencyDomain &fd) : Medium(fd, 1.0, 1.0)
    {
    }
    std::vector<real> layer_thicknesses(const std::vector<real> &in)
    {
        std::vector<real> ans(in.size() - 1);
        for (size_t i = 1; i < in.size(); i++)
        {
            assert(in[i] > in[i - 1]);
            ans[i - 1] = in[i] - in[i - 1];
        }

        return ans;
    }

    LayeredMedium::LayeredMedium(const FrequencyDomain &fd,
                                 const std::vector<real> &z_interfaces,
                                 const std::vector<Medium> &media,
                                 BC bottom_bc, BC top_bc)
        : fd(fd),
          z(z_interfaces),
          d(layer_thicknesses(z_interfaces)),
          media(media),
          bottom_bc(bottom_bc),
          top_bc(top_bc)
    {
    }

    int LayeredMedium::identify_layer(real z) const
    {
        assert(z < this->z.back());
        assert(z >= this->z.front());

        for (size_t n = 0; n < media.size(); n++)
            if (z >= this->z[n] && z < this->z[n + 1])
                return static_cast<int>(n);

        throw std::runtime_error("LayeredMedium::identify_layer failed.");
    }

    HalfSpace::HalfSpace(const FrequencyDomain &fd, const Medium &ground)
        : LayeredMedium(fd,
                        std::vector<real>{-INFINITY, 0, INFINITY},
                        std::vector<Medium>{ground, Vacuum(fd)},
                        BC::open,
                        BC::open)
    {
    }

    HalfSpacePEC::HalfSpacePEC(const FrequencyDomain &fd)
        : LayeredMedium(fd,
                        std::vector<real>{0, INFINITY},
                        std::vector<Medium>{Vacuum(fd)},
                        BC::PEC,
                        BC::open)
    {
    }

    HalfSpacePMC::HalfSpacePMC(const FrequencyDomain &fd)
        : LayeredMedium(fd,
                        std::vector<real>{0, INFINITY},
                        std::vector<Medium>{Vacuum(fd)},
                        BC::PMC,
                        BC::open)
    {
    }

} // namespace mthesis
