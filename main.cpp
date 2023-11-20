#include "fluxes.h"
#include "iostream"

int main(){
    double f = flux(fromm, 10.0, 9.0, 8.0, 7.0, 1.0, 0.1, 0.2);
    std::cout
        << f << std::endl;
}
