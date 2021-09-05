#ifndef MTHESIS_LEGENDRE_P_RECURRENCE_HPP
#define MTHESIS_LEGENDRE_P_RECURRENCE_HPP

#include <complex>
#include <iostream>

namespace mthesis {

template<typename T>
T legendre_p_recurrence(unsigned nu, T z)
{
    if (nu > 30) {
        std::cerr << "WARNING: legendre_p_recurrence "
                     "is slow for large degrees as nu = " << nu << "\n";
    }
    // See NIST Handbook of Mathematical Functions (14.10.3).
    T val;

    switch (nu) {
    case 0:
        val = 1;
        break;
    case 1:
        val = z;
        break;
    default:
        T t1 = static_cast<T>(2 * nu - 1) * z * legendre_p_recurrence(nu - 1, z);
        T t2 = static_cast<T>(nu - 1) * legendre_p_recurrence(nu - 2, z);
        T den = static_cast<T>(nu);
        val = (t1 - t2) / den;
        break;
    }

    return val;
}

} // namespace mthesis

#endif
