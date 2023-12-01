#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include "fluxes.h"
#include "pressures.h"
#include "grid.h"

typedef std::vector<std::vector<double>> Mat;

class idealSimulation : grid
{
public:
    idealSimulation(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY);
    void update(int nSteps);
    void xUpdate(int nSteps);
    using grid::print;

private:
    double gamma;

    Mat pressure();
    void xSources(Mat p);
    void ySources(Mat p);
    void updateDeltaT(Mat p);
};