#pragma once
#include <vector>
#include <iostream>
#include <iomanip>

typedef std::vector<std::vector<double>> Mat;

class grid
{
public:
    grid(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY);
    void advUpdate(int nSteps);
    void update(int nSteps);
    void xUpdate(int nSteps);
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
    void xAdvection(double func (double));
    void yAdvection(double func (double));
    Mat pressure();
    Mat uPressure();
    Mat vPressure();
    void xSources(Mat p);
    void ySources(Mat p);
    void xBorderCondition();
    void yBorderCondition();
    Mat pBorderCondition(Mat p);
};