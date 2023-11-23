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
    double deltaT;
    Mat Q1;
    Mat Q2x;
    Mat Q2y;
    Mat Q3;
    void xSweep();
    void ySweep();
    Mat pressure();
    void xSources(Mat p);
    void ySources(Mat p);
    void xBorderCondition();
    void yBorderCondition();
};