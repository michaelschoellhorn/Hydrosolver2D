#include "grid.h"
#include "fluxes.h"
#include "pressures.h"

grid::grid(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY)
{
    this->deltaX = deltaX;
    this->deltaY = deltaY;
    this->deltaT = 0.02;
    if (!QOne.empty())
    {
        // Get number of rows
        int sx = 0;
        for (const auto &row : QOne)
        {
            sx += 1;
        }

        // Get number of cols
        activYCells = int(QOne[0].size());

        activXCells = sx;
        ghostCells = 2;
        Nx = activXCells + 2 * ghostCells;
        Ny = activYCells + 2 * ghostCells;

        // Initialize grid matrizes
        Q1 = Mat(Nx, std::vector<double>(Ny, 0.0));
        Q2x = Mat(Nx, std::vector<double>(Ny, 0.0));
        Q2y = Mat(Nx, std::vector<double>(Ny, 0.0));
        Q3 = Mat(Nx, std::vector<double>(Ny, 0.0));

        // Simple check if dimensions are equal
        if ((activYCells != QTwox[0].size()) || (activYCells != QTwoy[0].size()) || (activYCells != QThree[0].size()))
        {
            std::cout << "Dimensions of QOne QTwo and QThree don't match!" << std::endl;
        }
        else
        {
            for (int i = ghostCells; i < activXCells + ghostCells; ++i)
            {
                for (int j = ghostCells; j < activYCells + ghostCells; ++j)
                {
                    // Set values to the grid
                    Q1[i][j] = QOne[i - ghostCells][j - ghostCells];
                    Q2x[i][j] = QTwox[i - ghostCells][j - ghostCells];
                    Q2y[i][j] = QTwoy[i - ghostCells][j - ghostCells];
                    Q3[i][j] = QThree[i - ghostCells][j - ghostCells];
                }
            }
        }
    }
    else
    {
        std::cout << "QOne is empty" << std::endl;
    }
}

void grid::print()
{
    // Q1
    std::cout << "Q1: \n";
    for (const auto &element : Q1)
    {
        for (const auto elem : element)
        {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    // Q2x
    std::cout << "Q2x: \n";
    for (const auto &element : Q2x)
    {
        for (const auto elem : element)
        {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

}

void grid::xBorderCondition()
{
    for (int x = 0; x < ghostCells; ++x)
    {
        for (int y = ghostCells; y < Ny - ghostCells; ++y)
        {
            Q1[x][y] = Q1[activXCells + x][y];
            Q2x[x][y] = Q2x[activXCells + x][y];
            Q2y[x][y] = Q2y[activXCells + x][y];
            Q3[x][y] = Q3[activXCells + x][y];

            Q1[Nx - ghostCells + x][y] = Q1[ghostCells + x][y];
            Q2x[Nx - ghostCells + x][y] = Q2x[ghostCells + x][y];
            Q2y[Nx - ghostCells + x][y] = Q2y[ghostCells + x][y];
            Q3[Nx - ghostCells + x][y] = Q3[ghostCells + x][y];
        }
    }
}

void grid::yBorderCondition()
{
    for (int y = 0; y < ghostCells; ++y)
    {
        for (int x = ghostCells; x < Nx - ghostCells; ++x)
        {
            Q1[x][y] = Q1[x][activYCells + y];
            Q2x[x][y] = Q2x[x][activYCells + y];
            Q2y[x][y] = Q2y[x][activYCells + y];
            Q3[x][y] = Q3[x][activYCells + y];

            Q1[x][Ny - ghostCells + y] = Q1[x][ghostCells + y];
            Q2x[x][Ny - ghostCells + y] = Q2x[x][ghostCells + y];
            Q2y[x][Ny - ghostCells + y] = Q2y[x][ghostCells + y];
            Q3[x][Ny - ghostCells + y] = Q3[x][ghostCells + y];
        }
    }
}

Mat grid::pBorderCondition(Mat p)
{
    for (int y = 0; y < ghostCells; ++y)
    {
        for (int x = ghostCells; x < Nx - ghostCells; ++x)
        {
            p[x][y] = p[x][activYCells + y];
            p[x][Ny - ghostCells + y] = p[x][ghostCells + y];
        }
    }
    for (int x = 0; x < ghostCells; ++x)
    {
        for (int y = ghostCells; y < Ny - ghostCells; ++y)
        {
            p[x][y] = p[activXCells + x][y];
            p[Nx - ghostCells + x][y] = p[ghostCells + x][y];
        }
    }
    return p;
}

void grid::update(int nSteps)
{
    xBorderCondition();
    yBorderCondition();
    for (size_t i = 0; i < nSteps; i++)
    {
        // xSweep ySweep
        xAdvection(minMod);
        Mat p = uPressure();
        xSources(p);
        yAdvection(minMod);
        p = vPressure();
        ySources(p);

        // ySweep xSweep
        yAdvection(minMod);
        p = vPressure();
        ySources(p);
        xAdvection(minMod);
        p = uPressure();
        xSources(p);
        print();
    }
}

void grid::advUpdate(int nSteps)
{
    xBorderCondition();
    yBorderCondition();
    for (size_t i = 0; i < nSteps; i++)
    {
        // xSweep ySweep
        xAdvection(minMod);
        yAdvection(minMod);

        // ySweep xSweep
        yAdvection(minMod);
        xAdvection(minMod);
        print();
    }
}

void grid::xAdvection(double func(double))
{
    // Initialize flux vectors
    std::vector<double> F1(Nx, 0.0);
    std::vector<double> F2x(Nx, 0.0);
    std::vector<double> F2y(Nx, 0.0);
    std::vector<double> F3(Nx, 0.0);

    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        // calculating fluxes
        double ux;
        for (int x = ghostCells; x < activXCells + ghostCells + 1; ++x)
        {
            ux = Q2x[x][y] / (Q1[x][y] + 1E-16);
            F1[x] = flux(func, Q1[x - 2][y], Q1[x - 1][y], Q1[x][y], Q1[x + 1][y], ux, deltaT, deltaX);
            F2x[x] = flux(func, Q2x[x - 2][y], Q2x[x - 1][y], Q2x[x][y], Q2x[x + 1][y], ux, deltaT, deltaX);
            F2y[x] = flux(func, Q2y[x - 2][y], Q2y[x - 1][y], Q2y[x][y], Q2y[x + 1][y], ux, deltaT, deltaX);
            F3[x] = flux(func, Q3[x - 2][y], Q3[x - 1][y], Q3[x][y], Q3[x + 1][y], ux, deltaT, deltaX);
        }
        // calculating Qhalf
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {

            Q1[x][y] += deltaT / deltaX * (F1[x] - F1[x + 1]);
            Q2x[x][y] += deltaT / deltaX * (F2x[x] - F2x[x + 1]);
            Q2y[x][y] += deltaT / deltaX * (F2y[x] - F2y[x + 1]);
            Q3[x][y] += deltaT / deltaX * (F3[x] - F3[x + 1]);
        }
    }
    xBorderCondition();
    yBorderCondition();
}

void grid::yAdvection(double func(double))
{
    // Initialize flux vectors
    std::vector<double> F1(Ny, 0.0);
    std::vector<double> F2x(Ny, 0.0);
    std::vector<double> F2y(Ny, 0.0);
    std::vector<double> F3(Ny, 0.0);

    for (int x = ghostCells; x < activXCells + ghostCells; ++x)
    {
        // calculating fluxes
        double uy;
        for (int y = ghostCells; y < activYCells + ghostCells + 1; ++y)
        {
            uy = Q2y[x][y] / (Q1[x][y] + 1E-16);
            F1[y] = flux(func, Q1[x][y - 2], Q1[x][y - 1], Q1[x][y], Q1[x][y + 1], uy, deltaT, deltaY);
            F2x[y] = flux(func, Q2x[x][y - 2], Q2x[x][y - 1], Q2x[x][y], Q2x[x][y + 1], uy, deltaT, deltaY);
            F2y[y] = flux(func, Q2y[x][y - 2], Q2y[x][y - 1], Q2y[x][y], Q2y[x][y + 1], uy, deltaT, deltaY);
            F3[y] = flux(func, Q3[x][y - 2], Q3[x][y - 1], Q3[x][y], Q3[x][y + 1], uy, deltaT, deltaY);
        }
        // calculating Qhalf
        for (int y = ghostCells; y < activYCells + ghostCells; ++y)
        {

            Q1[x][y] += deltaT / deltaY * (F1[y] - F1[y + 1]);
            Q2x[x][y] += deltaT / deltaY * (F2x[y] - F2x[y + 1]);
            Q2y[x][y] += deltaT / deltaY * (F2y[y] - F2y[y + 1]);
            Q3[x][y] += deltaT / deltaY * (F3[y] - F3[y + 1]);
        }
    }
    xBorderCondition();
    yBorderCondition();
}

Mat grid::pressure()
{
    Mat p(Nx, std::vector<double>(Ny, 0.0));
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            p[x][y] = idealPressure(1.4, Q1[x][y], Q2x[x][y], Q2y[x][y], Q3[x][y]);
        }
    }
    return pBorderCondition(p);
}

Mat grid::uPressure()
{
    Mat p(Nx, std::vector<double>(Ny, 0.0));
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            p[x][y] = viscIdealPressure(2.0, 1.4, Q1[x][y], Q2x[x][y], Q2y[x][y], Q3[x][y], Q1[x - 1][y], Q2x[x - 1][y], Q1[x + 1][y], Q2x[x + 1][y]);
        }
    }
    return pBorderCondition(p);
}

Mat grid::vPressure()
{
    Mat p(Nx, std::vector<double>(Ny, 0.0));
    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {
            p[x][y] = viscIdealPressure(2.0, 1.4, Q1[x][y], Q2x[x][y], Q2y[x][y], Q3[x][y], Q1[x][y - 1], Q2y[x][y - 1], Q1[x][y + 1], Q2y[x][y + 1]);
        }
    }
    return pBorderCondition(p);
}

void grid::xSources(Mat p)
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

void grid::ySources(Mat p)
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