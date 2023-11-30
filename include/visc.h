#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include "fluxes.h"
#include "pressures.h"
#include "grid.h"

class viscSimulation : public grid
{
public:
    viscSimulation(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY, double visc);
    void update(int nSteps);
    void xUpdate(int nSteps);

private:
    double gamma;
    double visc;

    Mat uPressure();
    Mat vPressure();
    void xSources(Mat p);
    void ySources(Mat p);
    void updateDeltaT(Mat px, Mat py);
};