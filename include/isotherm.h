#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include "fluxes.h"
#include "pressures.h"
#include "grid.h"

typedef std::vector<std::vector<double>> Mat;

class isothermalSimulation : public grid
{
public:
    isothermalSimulation(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY);
    void update(int nSteps);
    void xUpdate(int nSteps);

private:
    double cs;
    Mat pressure();
    void xSources(Mat p);
    void ySources(Mat p);
    void updateDeltaT();
};