#include "grid.h"

grid::grid(Mat QOne, Mat QTwo, Mat QThree, double deltaX, double deltaY)
{
    this->deltaX = deltaX;
    this->deltaY = deltaY;
    if (!QOne.empty())
    {
        // Get number of rows
        int sx = 0;
        for (const auto& row : QOne){
            sx += 1;
        }

        // Get number of cols
        activYCells = int(QOne[0].size());

        activXCells = sx;
        ghostCells = 2;

        // Initialize grid matrizes
        Q1 = Mat(activXCells + 2 * ghostCells, std::vector<double>(activYCells + 2 * ghostCells, 0.0));
        Q2 = Mat(activXCells + 2 * ghostCells, std::vector<double>(activYCells + 2 * ghostCells, 0.0));
        Q3 = Mat(activXCells + 2 * ghostCells, std::vector<double>(activYCells + 2 * ghostCells, 0.0));

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
                    Q1[i][j] = QOne[i-ghostCells][j-ghostCells];
                    Q2[i][j] = QTwo[i-ghostCells][j-ghostCells];
                    Q3[i][j] = QThree[i-ghostCells][j-ghostCells];
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
        for (const auto elem: element){
        std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

void grid::update() {}

void grid::advection1D() {}

void grid::sources1D() {}