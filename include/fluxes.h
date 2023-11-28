# pragma once
#include <cmath>
#include <algorithm>

inline double heavySide(const double x)
{
    if (x >= 0.0)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

inline double r(double Qmin2, double Qmin1, double Q, double Qp1, double u)
{
    return heavySide(u) * (Qmin1 - Qmin2) / (Q - Qmin1 + 1E-16) + heavySide(-u) * (Qp1 - Q) / (Q - Qmin1 + 1E-16);
}

inline double donorCell(double r)
{
    return 0.0;
}

inline double laxWendroff(double r)
{
    return 1.0;
}

inline double beamWarming(double r)
{
    return r;
}

inline double fromm(double r)
{
    return 0.5 * (1.0 + r);
}

inline double minMod(double r)
{
    if (r > 0.0)
    {
        if (r > 1.0)
        {
            return 1.0;
        }
        else
        {
            return r;
        }
    }
    else
    {
        return 0.0;
    }
}

inline double superBee(double r)
{
    return std::max({0.0, std::min(1.0, 2 * r), std::min(2.0, r)});
}

inline double MC(double r)
{
    return std::max(0.0, std::min({(1.0 + r) / 2.0, 2.0, 2.0 * r}));
}

inline double vanLeer(double r)
{
    return (r + std::abs(r)) / (1 + std::abs(r));
}

inline double flux(double function(double), double Qmin2, double Qmin1, double Q, double Qp1, double u, double deltaT, double deltaX)
{
    return 0.5 * u * (2.0 * heavySide(u) * Qmin1 + 2 * heavySide(-u) * Q) + 0.5 * std::abs(u) * (1.0 - std::abs(u * deltaT / deltaX)) * function(r(Qmin2, Qmin1, Q, Qp1, u)) * (Q - Qmin1);
}