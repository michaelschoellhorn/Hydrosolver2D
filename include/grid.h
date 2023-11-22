#pragma once
#include <vector>
#include <iostream>

typedef std::vector<std::vector<double>> Mat;

class grid
{
public:
    grid(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY);
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
    Mat Q2x;
    Mat Q2y;
    Mat Q3;
    void xSweep();
    void ySweep();
    void sources1D();
    void xBorderCondition();
    void yBorderCondition();
};