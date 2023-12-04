#include "fluxes.h"
#include "visc.h"
#include "ideal.h"
#include "isotherm.h"
#include "loadData.h"

int main()
{
    /*
    Mat QOne(40, std::vector<double>(40, 1.0));
    for (size_t i = 0; i < 20; i++)
    {
        for (size_t j = 0; j < 20; j++)
        {
            QOne[i][j] = 2.0;
        }
        
    }
    Mat QTwox(40, std::vector<double>(40, 0.0));
    Mat QTwoy(40, std::vector<double>(40, 0.0));
    Mat QThree(40, std::vector<double>(40, 1.0));
    for (size_t i = 0; i < 20; i++)
    {
        for (size_t j = 0; j < 20; j++)
        {
            QThree[i][j] = 2.0;
        }
    }
    */
    Mat QOne = loadFromTxt("startingRho.txt");
    Mat QTwox = loadFromTxt("startingRhoUx.txt");
    Mat QTwoy = loadFromTxt("startingRhoUy.txt");
    Mat QThree = loadFromTxt("startingRhoEps.txt");
    viscSimulation A(QOne, QTwox, QTwoy, QThree, 0.1, 0.1, 2.0);
    A.print();
    A.update(20);
    A.update(20);
}
