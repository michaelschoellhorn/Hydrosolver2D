#pragma once
#include <vector>
#include <iostream>

typedef std::vector<std::vector<double>> Mat;

class grid
{
public:
    grid(Mat QOne, Mat QTwo, Mat QThree, double deltaX, double deltaY);
    void update();
    void print();

private:
    int activXCells;
    int activYCells;
    int ghostCells;
    int Nx;
    int Ny;
    double deltaX;
    double deltaY;
    Mat Q1;
    Mat Q2;
    Mat Q3;
    void advection1D();
    void sources1D();
    void borderCondition();
};