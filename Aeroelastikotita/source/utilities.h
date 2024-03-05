#pragma once

//#include <Eigen/src/Core/util/Constants.h>
#include <complex>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <string>
#include <numbers>
#include <vector>
#include <iomanip>

#define dt 0.005

const int parInitVal = 25;
const double pi = std::numbers::pi;

struct cmpx {
    double len, theta, real, imag;
    void set(std::complex<double> val) {
        len = sqrt(pow(val.real(),2) + pow(val.imag(),2));
        theta = atan2(val.imag(), val.real());
        real = val.real();
        imag = val.imag();
    }
};

struct inputData {
    // Struct variables
    int siz;
    double alpha, theta, U_inf, rho, chord, A1, A2, B1, B2, t_tot;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M, K, C;
    Eigen::Matrix<double, Eigen::Dynamic, 2> Lift, Drag;
    Eigen::Matrix<cmpx,4,2> eigenvectors;
    Eigen::Matrix<cmpx,4,2> prevEigenvectors;
    Eigen::Vector<double,2> Q;
    Eigen::Vector<double,2> prevQ;
    Eigen::Vector<double,2> partialSol;

    std::vector<std::complex<double> > eigenvalues;
    std::vector<std::complex<double> > prevEigenvalues;

    // Methods declaration
    double getCD(double a);
    double getCDgrad(double a);
    double getCL(double a);
    double getCLgrad(double a);
    void readInput(std::string filename);
    void initialize();
    void calcDamp();
    void computeEigen();
    void computeConstants();
    void calcNextStep();
    void calcSteadyLoads();
    void nextAnalyticalStep(std::string filename);
    double getAeff();


    // State variables
    double u, uv, ua, w, wv, wa;

    bool theodor; // Unsteady solution flag (True)
    
    inputData(int sz)
    {
        siz = sz;
        M.resize(sz,sz);
        K.resize(sz,sz);
        C.resize(sz,sz);
        M.setZero(); K.setZero(); C.setZero();
        Lift.resize(parInitVal,2);
        Drag.resize(parInitVal,2);
        eigenvalues.resize(2);
        prevEigenvalues.resize(2);
        prevQ.setZero();
        u = 0;
        uv = 0;
        ua = 0;
        w = 0;
        wv = 0;
        wa = 0;
        t_tot = 0;
        Q.setZero();
        theodor = false;
    }
};

// Struct used to compute unsteady variants 
struct Theodorsen {
    inputData *target;
    double b[2];
    double A[2];
    double y[2];
    double ae, L, D;
    double grad0, a0;

    void calcNextY();
    void calcDamp();
    void calcLoads();
    // Constructor assigns target struct pointer
    Theodorsen(inputData* p) : target(p){
        b[0] = 0.0455;
        b[1] = 0.3;
        A[0] = 0.165;
        A[1] = 0.335;
        if (target->t_tot != 0) {
            std::cout << "Instantiated Theodorsen on non-initial solution" << std::endl;
        }else {
            y[0] = target->getAeff()*A[0];
            y[1] = target->getAeff()*A[1];
        }

        grad0 = (target->getCL(0.0+5.0*pi/180.0)-target->getCL(0.0-5.0*pi/180.0))/((double)2.0*5.0/180.0*pi);

        if (grad0 != 0)
        {
            a0 = -target->getCL(0)/grad0;
        }else {
            std::cout << "zero CL gradient" << std::endl;
            a0 = 0;
        }
    };
};

void readAero(std::string filenameCl, std::string filenameCd, inputData &input);

void readMat(int sz, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &M, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &K, std::string filename);
