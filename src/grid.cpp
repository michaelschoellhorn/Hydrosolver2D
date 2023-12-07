#include "grid.h"

grid::grid(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY)
{
    this->deltaX = deltaX;
    this->deltaY = deltaY;
    cfl = 0.5;
    this->deltaT = 0.005;
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

void grid::saveTo(std::string saveFile)
{
    std::ofstream outFile(saveFile);
    if (!outFile.is_open()) // Check if file is open
    {
        std::cerr << "Can't open saveFile" << saveFile << std::endl; // Error
    }
    outFile << "NumX Numy Ghostcells CFL deltaX deltaY\n";
    outFile << activXCells << " " << activYCells << " " << ghostCells << " " << cfl << " " << deltaX << " " << deltaY << "\n" << std::endl;
    // Print Rho
    for (const auto &row : Q1)
    {
        bool isFirst = 1;
        for (const auto &value : row)
        {
            if (!isFirst)
            {
                outFile << " ";
            }
            outFile << value;
            isFirst = 0;
        }
        outFile << "\n";
    }
    outFile << std::endl;
    // Print Rho ux
    for (const auto &row : Q2x)
    {
        bool isFirst = 1;
        for (const auto &value : row)
        {
            if (!isFirst)
            {
                outFile << " ";
            }
            outFile << value;
            isFirst = 0;
        }
        outFile << "\n";
    }
    outFile << std::endl;
    // Print Rho uy
    for (const auto &row : Q2y)
    {
        bool isFirst = 1;
        for (const auto &value : row)
        {
            if (!isFirst)
            {
                outFile << " ";
            }
            outFile << value;
            isFirst = 0;
        }
        outFile << "\n";
    }
    outFile << std::endl;
    // Print Rho epsilon
    for (const auto &row : Q3)
    {
        bool isFirst = 1;
        for (const auto &value : row)
        {
            if (!isFirst)
            {
                outFile << " ";
            }
            outFile << value;
            isFirst = 0;
        }
        outFile << "\n";
    }
    outFile << std::endl;
}

void grid::print()
{
    // Q1
    std::cout << "Q1: \n";
    for (const auto &element : Q1)
    {
        for (const auto elem : element)
        {
            std::cout << std::fixed << std::setprecision(2) << elem << " ";
        }
        std::cout << std::endl;
    }
    // Q2x
    std::cout << "Q2x: \n";
    for (const auto &element : Q2x)
    {
        for (const auto elem : element)
        {
            std::cout << std::fixed << std::setprecision(2) << elem << " ";
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
        //print();
    }
}

void grid::xAdvection(double func(double))
{
    // Initialize flux vectors
    Mat F1(Nx, std::vector<double>(Ny, 0.0));
    Mat F2x(Nx, std::vector<double>(Ny, 0.0));
    Mat F2y(Nx, std::vector<double>(Ny, 0.0));
    Mat F3(Nx, std::vector<double>(Ny, 0.0));

    for (int y = ghostCells; y < activYCells + ghostCells; ++y)
    {
        // calculating fluxes
        double ux;
        for (int x = ghostCells; x < activXCells + ghostCells + 1; ++x)
        {
            ux = 0.5 * (Q2x[x][y] / (Q1[x][y] + 1E-16) + Q2x[x - 1][y] / (Q1[x - 1][y] + 1E-16));
            F1[x][y] = flux(func, Q1[x - 2][y], Q1[x - 1][y], Q1[x][y], Q1[x + 1][y], ux, deltaT, deltaX);
            F2x[x][y] = flux(func, Q2x[x - 2][y], Q2x[x - 1][y], Q2x[x][y], Q2x[x + 1][y], ux, deltaT, deltaX);
            F2y[x][y] = flux(func, Q2y[x - 2][y], Q2y[x - 1][y], Q2y[x][y], Q2y[x + 1][y], ux, deltaT, deltaX);
            F3[x][y] = flux(func, Q3[x - 2][y], Q3[x - 1][y], Q3[x][y], Q3[x + 1][y], ux, deltaT, deltaX);
        }
        // calculating Qhalf
        for (int x = ghostCells; x < activXCells + ghostCells; ++x)
        {

            Q1[x][y] += deltaT / deltaX * (F1[x][y] - F1[x + 1][y]);
            Q2x[x][y] += deltaT / deltaX * (F2x[x][y] - F2x[x + 1][y]);
            Q2y[x][y] += deltaT / deltaX * (F2y[x][y] - F2y[x + 1][y]);
            Q3[x][y] += deltaT / deltaX * (F3[x][y] - F3[x + 1][y]);
        }
    }
    xBorderCondition();
    yBorderCondition();
}

void grid::yAdvection(double func(double))
{
    // Initialize flux vectors
    Mat F1(Nx, std::vector<double>(Ny, 0.0));
    Mat F2x(Nx, std::vector<double>(Ny, 0.0));
    Mat F2y(Nx, std::vector<double>(Ny, 0.0));
    Mat F3(Nx, std::vector<double>(Ny, 0.0));

    for (int x = ghostCells; x < activXCells + ghostCells; ++x)
    {
        // calculating fluxes
        double uy;
        for (int y = ghostCells; y < activYCells + ghostCells + 1; ++y)
        {
            uy = 0.5 * (Q2y[x][y] / (Q1[x][y] + 1E-16) + Q2y[x][y - 1] / (Q1[x][y - 1] + 1E-16));
            F1[x][y] = flux(func, Q1[x][y - 2], Q1[x][y - 1], Q1[x][y], Q1[x][y + 1], uy, deltaT, deltaY);
            F2x[x][y] = flux(func, Q2x[x][y - 2], Q2x[x][y - 1], Q2x[x][y], Q2x[x][y + 1], uy, deltaT, deltaY);
            F2y[x][y] = flux(func, Q2y[x][y - 2], Q2y[x][y - 1], Q2y[x][y], Q2y[x][y + 1], uy, deltaT, deltaY);
            F3[x][y] = flux(func, Q3[x][y - 2], Q3[x][y - 1], Q3[x][y], Q3[x][y + 1], uy, deltaT, deltaY);
        }
        // calculating Qhalf
        for (int y = ghostCells; y < activYCells + ghostCells; ++y)
        {

            Q1[x][y] += deltaT / deltaY * (F1[x][y] - F1[x][y + 1]);
            Q2x[x][y] += deltaT / deltaY * (F2x[x][y] - F2x[x][y + 1]);
            Q2y[x][y] += deltaT / deltaY * (F2y[x][y] - F2y[x][y + 1]);
            Q3[x][y] += deltaT / deltaY * (F3[x][y] - F3[x][y + 1]);
        }
    }
    xBorderCondition();
    yBorderCondition();
}
