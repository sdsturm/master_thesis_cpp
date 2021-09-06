#ifndef MTHESIS_LEGENDRE_P_RECURRENCE_HPP
#define MTHESIS_LEGENDRE_P_RECURRENCE_HPP

#include <vector>

namespace mthesis {

template<typename  T>
T legendre_p_recurrence(unsigned nu, T z)
{
    // See (14.10.3) in NIST Handbook of Mathematical Functions.
    std::vector<T> P(nu + 1);

    P[0] = static_cast<T>(1);
    if (nu > 0) {
        P[1] = z;
        if (nu > 1) {
            for (unsigned l = 2; l <= nu; l++) {
                T den = static_cast<T>(l);
                T t1 = static_cast<T>(2 * l - 1) * z * P[l - 1];
                T t2 = static_cast<T>(l - 1) * P[l - 2];
                P[l] = (t1 - t2) / den;
            }
        }
    }
    return P.back();
}

} // namespace mthesis

#endif
