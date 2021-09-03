#include <mthesis/solution_domain.hpp>

#include <gsl/gsl_const_mksa.h>

#include <cassert>

namespace mthesis {

FrequencyDomain::FrequencyDomain(real f)
    : f(f),
      omega(2.0 * M_PI * f),
      k_0(omega / GSL_CONST_MKSA_SPEED_OF_LIGHT),
      lambda_0(2.0 * M_PI / k_0)
{
    assert(f > 0);
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
    assert(std::isfinite(std::abs(eps_r)));
    assert(std::isfinite(std::abs(mu_r)));
    assert(eps_r.real() > 0.0);
    assert(mu_r.imag() <= 0.0);
    assert(mu_r.real() > 0.0);
    assert(mu_r.imag() <= 0.0);
}

cmplx cmplx_permittivity(const FrequencyDomain &fd, real eps_r, real sigma)
{
    return cmplx(eps_r,
                 -sigma / (fd.omega * GSL_CONST_MKSA_VACUUM_PERMITTIVITY));
}

Vacuum::Vacuum(const FrequencyDomain &fd) : Medium(fd, 1.0, 1.0)
{
}

std::vector<real> layer_thicknesses(const std::vector<real> &z_interfaces)
{
    std::vector<real> d(z_interfaces.size() - 1);
    for (size_t i = 1; i < z_interfaces.size(); i++)
    {
        assert(z_interfaces[i] > z_interfaces[i - 1]);
        d[i - 1] = z_interfaces[i] - z_interfaces[i - 1];
    }
    return d;
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
    assert(z_interfaces.size() >= 2);
    assert(z_interfaces.size() == media.size() + 1);

    if (top_bc == BC::open)
        assert(std::isinf(z.back()));

    if (bottom_bc == BC::open)
        assert(std::isinf(z.front()));

    for (size_t n = 1; n < d.size() - 1; n++)
        assert(std::isfinite(d[n]));
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

FreeSpace::FreeSpace(const FrequencyDomain &fd)
    : LayeredMedium(fd,
                    std::vector<real>{-INFINITY, INFINITY},
                    std::vector<Medium>{Vacuum(fd)},
                    BC::open,
                    BC::open)
{
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
