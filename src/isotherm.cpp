#include "isotherm.h"

isothermalSimulation::isothermalSimulation(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY) : grid(QOne, QTwox, QTwoy, QThree, deltaX, deltaY)
{
    cs = 1.0;
}

void isothermalSimulation::update(int nSteps)
{
    xBorderCondition();
    yBorderCondition();
    for (size_t i = 0; i < nSteps; i++)
    {
        // xSweep ySweep
        xAdvection(minMod);
        Mat p = pressure();
        xSources(p);
        yAdvection(minMod);
        p = pressure();
        ySources(p);

        // ySweep xSweep
        yAdvection(minMod);
        p = pressure();
        ySources(p);
        xAdvection(minMod);
        p = pressure();
        xSources(p);
        updateDeltaT();
        // print();
    }
    print();
}

void isothermalSimulation::xUpdate(int nSteps)
{
    xBorderCondition();
    yBorderCondition();
    for (size_t i = 0; i < nSteps; i++)
    {
        // xSweep
        xAdvection(minMod);
        Mat p = pressure();
        xSources(p);
    }
    print();
}

Mat isothermalSimulation::pressure()
{
    Mat p(Nx, std::vector<double>(Ny, 0.0));
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            p[x][y] = isothermPressure(1.0, Q1[x][y]);
        }
    }
    return pBorderCondition(p);
}

void isothermalSimulation::xSources(Mat p)
{
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            Q2x[x][y] -= deltaT / (2 * deltaX) * (p[x + 1][y] - p[x - 1][y]);
            Q3[x][y] -= deltaT / (2 * deltaX) * (p[x + 1][y] * Q2x[x + 1][y] / (Q1[x + 1][y] + 1E-16) - p[x - 1][y] * Q2x[x - 1][y] / (Q1[x - 1][y] + 1E-16));
        }
    }
    xBorderCondition();
    yBorderCondition();
}

void isothermalSimulation::ySources(Mat p)
{
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            Q2y[x][y] -= deltaT / (2 * deltaY) * (p[x][y + 1] - p[x][y - 1]);
            Q3[x][y] -= deltaT / (2 * deltaY) * (p[x][y + 1] * Q2y[x][y + 1] / (Q1[x][y + 1] + 1E-16) - p[x][y - 1] * Q2y[x][y - 1] / (Q1[x][y - 1] + 1E-16));
        }
    }
    xBorderCondition();
    yBorderCondition();
}

void isothermalSimulation::updateDeltaT()
{
    double nuX = 0.0;
    double nuY = 0.0;
    double temp1;
    double temp2;
    double temp3;
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
        {
            for (int x = ghostCells; x < activXCells + ghostCells; ++x)
            {
                temp1 = cs + abs(Q2x[x][y] / (Q1[x][y] + 1E-16)); // |cs_x| + |u_x|
                temp2 = cs + abs(Q2y[x][y] / (Q1[x][y] + 1E-16)); // |cs_y| + |u_y|
                if (nuX < temp1)
                {
                    nuX = temp1;
                }
                if (nuY < temp2)
                {
                    nuY = temp2;
                }
            }
        }
    deltaT = cfl * std::min(deltaX / nuX, deltaY / nuY);
}