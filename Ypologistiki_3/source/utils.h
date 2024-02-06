#pragma once

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <math.h>
#include <Eigen/Core>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <math.h>

#define pi 3.14159265359

void writeResults(int iter);
void discretize();
void calcS();
void initialize();
void calcCons();
void calcPrim();
void calcFlux();
Eigen::Matrix<double, 3, 3> asmLeft(double U, double dens, double C);
Eigen::Matrix<double, 3, 3> asmRight(double U, double dens, double C);
void roeFlux();
void updateBC();

