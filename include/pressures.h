#pragma once
#include <cmath>
#include <algorithm>

inline double idealPressure(double gamma, double rho, double rhoU, double rhoV, double rhoE)
{
    return (gamma - 1.0) * (rhoE - 0.5 * (pow(rhoU, 2.0) + pow(rhoV, 2.0)) / (rho + 1E-16));
}

inline double isothermPressure(double cs, double rho)
{
    return pow(cs, 2.0) * rho;
}

inline double viscIdealPressure(double visc, double gamma, double rho, double rhoU, double rhoV, double rhoE, double rhoMin1, double rhoUMin1, double rhoP1, double rhoUP1){
    double p = idealPressure(gamma, rho, rhoU, rhoV, rhoE);
    double uP1 = rhoUP1 / (rhoP1 + 1E-16);
    double uMin1 = rhoUMin1 / (rhoMin1 + 1E-16);
    if (abs(uP1) <= abs(uMin1))
    {
        p += 0.25 * pow(visc, 2.0) * pow(abs(uP1) - abs(uMin1), 2.0) * rho;
    }
    return p;
}
/*
// old pressure function, depreciated
inline double viscIdealPressure(double visc, double gamma, double rho, double rhoU, double rhoV, double rhoE, double rhoMin1, double rhoUMin1, double rhoVMin1, double rhoP1, double rhoUP1, double rhoVP1){
    double p = idealPressure(gamma, rho, rhoU, rhoV, rhoE);
    double uP1 = rhoUP1 / (rhoP1 + 1E-16);
    double vP1 = rhoVP1 / (rhoP1 + 1E-16);
    double uMin1 = rhoUMin1 / (rhoMin1 + 1E-16);
    double vMin1 = rhoVMin1 / (rhoMin1 + 1E-16);
    if ((pow(pow(uP1, 2.0) + pow(vP1, 2.0), 0.5)) <= (pow(pow(uMin1, 2.0) + pow(vMin1, 2.0), 0.5)))
    {
        p += 0.25 * pow(visc, 2.0) * pow(pow(pow(uP1, 2.0) + pow(vP1, 2.0), 0.5) - pow(pow(uMin1, 2.0) + pow(vMin1, 2.0), 0.5), 2.0);
    }
    return p;
}
*/