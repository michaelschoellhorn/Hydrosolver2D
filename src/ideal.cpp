#include "ideal.h"

idealSimulation::idealSimulation(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY) : grid(QOne, QTwox, QTwoy, QThree, deltaX, deltaY)
{
    gamma = 1.4;
}

void idealSimulation::update(int nSteps)
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
        //updateDeltaT(p);
        // print();
    }
    print();
}

void idealSimulation::xUpdate(int nSteps)
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

Mat idealSimulation::pressure()
{
    Mat p(Nx, std::vector<double>(Ny, 0.0));
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            p[x][y] = idealPressure(gamma, Q1[x][y], Q2x[x][y], Q2y[x][y], Q3[x][y]);
        }
    }
    return pBorderCondition(p);
}

void idealSimulation::xSources(Mat p)
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

void idealSimulation::ySources(Mat p)
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

void idealSimulation::updateDeltaT(Mat p)
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
                temp1 = std::sqrt(gamma * p[x][y] / (Q1[x][y] + 1E-16));
                temp2 = temp1 + abs(Q2x[x][y] / (Q1[x][y] + 1E-16)); // |cs| + |u_x|
                temp3 = temp1 + abs(Q2y[x][y] / (Q1[x][y] + 1E-16)); // |cs| + |u_y|;
                if (nuX < temp2)
                {
                    nuX = temp2;
                }
                if (nuY < temp3)
                {
                    nuY = temp3;
                }
            }
        }
    deltaT = cfl * std::min(deltaX / nuX, deltaY / nuY);
}