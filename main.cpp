#include "fluxes.h"
#include "visc.h"
#include "ideal.h"
#include "isotherm.h"

int main()
{
    double f = flux(fromm, 10.0, 9.0, 8.0, 7.0, 1.0, 0.1, 0.2);
    std::cout
        << f << std::endl;
    // Mat QOne({std::vector<double>({1.5, 2.0, 3.5, 4.0}), std::vector<double>({1.0, 2.0, 3.0, 4.0}), std::vector<double>({5.0, 6.0, 7.0, 8.0}), std::vector<double>({5.0, 6.0, 7.0, 8.0})});
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
    viscSimulation A(QOne, QTwox, QTwoy, QThree, 0.1, 0.1, 2.0);
    A.print();
    A.update(20);
    A.update(20);
}
