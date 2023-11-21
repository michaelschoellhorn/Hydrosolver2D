#include "grid.h"

grid::grid(Mat QOne, Mat QTwo, Mat QThree, double deltaX, double deltaY)
{
    this->deltaX = deltaX;
    this->deltaY = deltaY;
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
        Q2 = Mat(Nx, std::vector<double>(Ny, 0.0));
        Q3 = Mat(Nx, std::vector<double>(Ny, 0.0));

        // Simple check if dimensions are equal
        if ((activYCells != QTwo[0].size()) || (activYCells != QThree[0].size()))
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
                    Q2[i][j] = QTwo[i - ghostCells][j - ghostCells];
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
    // Function to print Q1
    std::cout << "Q1: \n";
    for (const auto &element : Q1)
    {
        for (const auto elem : element)
        {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

void grid::borderCondition()
{
    for (int x = 0; x < ghostCells; ++x)
    {
        for (int y = ghostCells; y < Ny - ghostCells; ++y)
        {
            Q1[x][y] = Q1[activXCells + x][y];
            Q2[x][y] = Q2[activXCells + x][y];
            Q3[x][y] = Q3[activXCells + x][y];

            Q1[Nx - ghostCells + x][y] = Q1[ghostCells + x][y];
            Q2[Nx - ghostCells + x][y] = Q2[ghostCells + x][y];
            Q3[Nx - ghostCells + x][y] = Q3[ghostCells + x][y];
        }
    }
    for (int y = 0; y < ghostCells; ++y)
    {
        for (int x = ghostCells; x < Nx - ghostCells; ++x)
        {
            Q1[x][y] = Q1[x][activYCells + y];
            Q2[x][y] = Q2[x][activYCells + y];
            Q3[x][y] = Q3[x][activYCells + y];

            Q1[x][Ny - ghostCells + y] = Q1[x][ghostCells + y];
            Q2[x][Ny - ghostCells + y] = Q2[x][ghostCells + y];
            Q3[x][Ny - ghostCells + y] = Q3[x][ghostCells + y];
        }
    }
}

void grid::update()
{
    borderCondition();
}

void grid::advection1D()
{
}

void grid::sources1D()
{
}