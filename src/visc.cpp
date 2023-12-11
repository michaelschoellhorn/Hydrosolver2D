#include "visc.h"

viscSimulation::viscSimulation(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY, double visc) : grid(QOne, QTwox, QTwoy, QThree, deltaX, deltaY)
{
    gamma = 1.4;
    this->visc = visc;
}

void viscSimulation::update(int nSteps)
{
    xBorderCondition();
    yBorderCondition();
    for (size_t i = 0; i < nSteps; i++)
    {
        // xSweep ySweep
        xAdvection(donorCell);
        Mat px = uPressure();
        xSources(px);
        yAdvection(donorCell);
        Mat py = vPressure();
        ySources(py);

        // ySweep xSweep
        yAdvection(donorCell);
        px = vPressure();
        ySources(px);
        xAdvection(donorCell);
        py = uPressure();
        xSources(py);
        // updateDeltaT(px, py);
        // print();
    }
    //print();
}

void viscSimulation::xUpdate(int nSteps)
{
    xBorderCondition();
    yBorderCondition();
    for (size_t i = 0; i < nSteps; i++)
    {
        // xSweep
        xAdvection(minMod);
        Mat p = uPressure();
        xSources(p);
    }
    //print();
}

Mat viscSimulation::uPressure()
{
    Mat p(Nx, std::vector<double>(Ny, 0.0));
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            p[x][y] = viscIdealPressure(visc, gamma, Q1[x][y], Q2x[x][y], Q2y[x][y], Q3[x][y], Q1[x - 1][y], Q2x[x - 1][y], Q1[x + 1][y], Q2x[x + 1][y]);
        }
    }
    return pBorderCondition(p);
}

Mat viscSimulation::vPressure()
{
    Mat p(Nx, std::vector<double>(Ny, 0.0));
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            p[x][y] = viscIdealPressure(visc, gamma, Q1[x][y], Q2x[x][y], Q2y[x][y], Q3[x][y], Q1[x][y - 1], Q2y[x][y - 1], Q1[x][y + 1], Q2y[x][y + 1]);
        }
    }
    return pBorderCondition(p);
}

void viscSimulation::xSources(Mat p)
{
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            Q3[x][y] -= deltaT / (2 * deltaX) * (p[x + 1][y] * Q2x[x + 1][y] / (Q1[x + 1][y] + 1E-16) - p[x - 1][y] * Q2x[x - 1][y] / (Q1[x - 1][y] + 1E-16));
        }
    }
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
        Q2x[x][y] -= deltaT / (2 * deltaX) * (p[x + 1][y] - p[x - 1][y]);
        }
    }
    xBorderCondition();
    yBorderCondition();
}

void viscSimulation::ySources(Mat p)
{
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            Q2y[x][y] -= deltaT / (2 * deltaY) * (p[x][y + 1] - p[x][y - 1]);
            Q3[x][y] -= deltaT / (2 * deltaY) * (p[x][y + 1] * Q2y[x][y + 1] / (Q1[x][y + 1] + 1E-16) - p[x][y - 1] * Q2y[x][y - 1] / (Q1[x][y - 1] + 1E-16));
        }
    }
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            Q2y[x][y] -= deltaT / (2 * deltaY) * (p[x][y + 1] - p[x][y - 1]);
            }
    }
    xBorderCondition();
    yBorderCondition();
}

void viscSimulation::updateDeltaT(Mat px, Mat py)
{
    double nuX = 0.0;
    double nuY = 0.0;
    double temp1;
    double temp2;
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            temp1 = std::sqrt(gamma * px[x][y] / (Q1[x][y] + 1E-16)) + abs(Q2x[x][y] / (Q1[x][y] + 1E-16)); // |cs_x| + |u_x|
            temp2 = std::sqrt(gamma * py[x][y] / (Q1[x][y] + 1E-16)) + abs(Q2y[x][y] / (Q1[x][y] + 1E-16)); // |cs_y| + |u_y|
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