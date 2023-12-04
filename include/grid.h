#pragma once
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "fluxes.h"
#include "pressures.h"

typedef std::vector<std::vector<double>> Mat;

class grid
{
public:
    grid(Mat QOne, Mat QTwox, Mat QTwoy, Mat QThree, double deltaX, double deltaY);
    void saveTo(std::string saveFile);
    void advUpdate(int nSteps);
    void print();

protected:
    int activXCells;
    int activYCells;
    int ghostCells;
    int Nx;
    int Ny;

    double cfl;
    double deltaX;
    double deltaY;
    double deltaT;

    Mat Q1;
    Mat Q2x;
    Mat Q2y;
    Mat Q3;

    void xAdvection(double func (double));
    void yAdvection(double func (double));

    void xBorderCondition();
    void yBorderCondition();
    Mat pBorderCondition(Mat p);
};