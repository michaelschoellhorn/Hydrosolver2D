#include "fluxes.h"
#include "grid.h"


int main(){
    double f = flux(fromm, 10.0, 9.0, 8.0, 7.0, 1.0, 0.1, 0.2);
    std::cout
        << f << std::endl;
    Mat QOne(4, std::vector<double>(4, 2.0));
    Mat QTwo(4, std::vector<double>(4, 0.0));
    Mat QThree(4, std::vector<double>(4, 0.0));
    grid A(QOne, QTwo, QThree, 0.1, 0.1);
    A.print();
}
