#include "fluxes.h"
#include "grid.h"

int main()
{
    double f = flux(fromm, 10.0, 9.0, 8.0, 7.0, 1.0, 0.1, 0.2);
    std::cout
        << f << std::endl;
    Mat QOne({std::vector<double>({-1.5, -2.0, -3.5, -4.0}), std::vector<double>({1.0, 2.0, 3.0, 4.0}), std::vector<double>({5.0, 6.0, 7.0, 8.0}), std::vector<double>({-5.0, -6.0, -7.0, -8.0})});
    Mat QTwox(4, std::vector<double>(4, 1.0));
    Mat QTwoy(4, std::vector<double>(4, 0.0));
    Mat QThree(4, std::vector<double>(4, 0.0));
    grid A(QOne, QTwox, QTwoy, QThree, 0.1, 0.1);
    A.print();
    A.update();
    A.print();
}
