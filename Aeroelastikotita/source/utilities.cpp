#include "utilities.h"
#include <Eigen/src/Eigenvalues/GeneralizedEigenSolver.h>
#include <fstream>
#include <istream>


// Reads mass and stiffness matrices
void readMat(int sz, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &M, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &K, std::string filename)
{
    std::ifstream mat(filename); 
    std::string a,b,c,trash;

    mat >> trash;

    for (int i = 0; i < sz; i++) {
        mat >> a >> b ;
        M(i,0) = std::stod(a);
        M(i,1) = std::stod(b);
    }

    mat >> trash;

    for (int i = 0; i < sz; i++) {
        mat >> a >> b ;
        K(i,0) = std::stod(a);
        K(i,1) = std::stod(b);
    }

}

// Reads input file
void inputData::readInput(std::string filename) {
    std::ifstream inp(filename);

    std::string trash, text;
    
    getline(inp,trash);
    getline(inp,text);
    U_inf = std::stod(text);

    getline(inp,trash);
    getline(inp,text);
    theta = std::stod(text)*pi/180;

    getline(inp,trash);
    getline(inp,text);
    alpha = std::stod(text)*pi/180;

    getline(inp,trash);
    getline(inp,text);
    rho = std::stod(text);

    getline(inp,trash);
    getline(inp,text);
    chord = std::stod(text);

    getline(inp,trash);
    getline(inp,text);
    u = std::stod(text);

    getline(inp,trash);
    getline(inp,text);
    uv = std::stod(text);

    getline(inp,trash);
    getline(inp,text);
    ua = std::stod(text);

    getline(inp,trash);
    getline(inp,text);
    w = std::stod(text);

    getline(inp,trash);
    getline(inp,text);
    wv = std::stod(text);

    getline(inp,trash);
    getline(inp,text);
    wa = std::stod(text);

    // MISPLACED UTILITY
    
    // Transform stiffness matrix
    Eigen::Matrix<double,2 ,2> R;
    R(0,0) = cos(theta);
    R(0,1) = sin(theta);
    R(1,0) = -sin(theta);
    R(1,1) = cos(theta);

    K = R.transpose()*K*R;
}

// Reads aerodynamic data (curves)
void readAero(std::string filenameCl, std::string filenameCd, inputData &input) {

    int Lsz;
    int Dsz;
    std::string textCl;
    std::string textCd;

    std::ifstream cl(filenameCl);
    std::ifstream cd(filenameCd);

    // First entry is the count of the data points
    cl >> Lsz;
    cd >> Dsz;

    input.Lift.conservativeResize(Lsz,2);
    input.Drag.conservativeResize(Dsz,2);

    for (int i = 0; i < Lsz; i++)
    {
        cl >> input.Lift(i,0) >> input.Lift(i,1);
    }

    for (int i = 0; i < Dsz; i++)
    {
        cd >> input.Drag(i,0) >> input.Drag(i,1);
    }
}

double inputData::getCL(double a /*Radians*/){
    double alph = a/pi*180;
    double cl;
    int i = 0; 
    int sz = Lift.rows()-1;

    if (alph > 50) alph = 50;
    if (alph < -50) alph = -50;

    // Find section for linear interpolation
    while (i < Lift.rows()) {
        if (alph <= Lift(i,0))
        {
            break;    
        }
        i++;
    }

    // If a is outside of data bounds
    // 2nd order extrapolation

    double grad0 = (Lift(1,1)-Lift(0,1))/(Lift(1,0)-Lift(0,0));
    double grad0_1 = (Lift(2,1)-Lift(0,1))/(Lift(2,0)-Lift(0,0));

    double gradn = (Lift(sz,1)-Lift(sz-1,1))/(Lift(sz,0)-Lift(sz-1,0));
    double gradn_1 = (Lift(sz,1)-Lift(sz-2,1))/(Lift(sz,0)-Lift(sz-2,0));

    double curv0 = (grad0_1-grad0)/(Lift(1,0)-Lift(0,0));
    double curvn = (gradn-gradn_1)/(Lift(sz,0)-Lift(sz-1,0));

    if (i==0) {
        // Extrapolation using curvature
        cl = Lift(0,1) + grad0*(alph-Lift(0,0)) + pow(alph-Lift(0,0),2)/2*curv0; 
    } else if (i == Lift.rows()) {
        // Extrapolation using curvature
        cl = Lift(sz,1) + gradn*(alph-Lift(sz,0)) + pow(alph-Lift(sz,0),2)/2*curvn; 
    }else {
        // linear interpolation
        cl = Lift(i-1,1) + (Lift(i,1)-Lift(i-1,1))/(Lift(i,0)-Lift(i-1,0))*(alph-Lift(i-1,0));
    }

    return cl;
}

double inputData::getCD(double a /*Radians*/){
    double alph = a/pi*180;
    double cd;
    int i = 0; 
    int sz = Drag.rows()-1;

    if (alph > 50) alph = 50;
    if (alph < -50) alph = -50;

    while (i < Drag.rows()) {
        if (alph <= Drag(i,0))
        {
            break;    
        }
        i++;
    }

    // If a is outside of data bounds
    // 2nd order extrapolation

    double grad0 = (Drag(1,1)-Drag(0,1))/(Drag(1,0)-Drag(0,0));
    double grad0_1 = (Drag(2,1)-Drag(0,1))/(Drag(2,0)-Drag(0,0));

    double gradn = (Drag(sz,1)-Drag(sz-1,1))/(Drag(sz,0)-Drag(sz-1,0));
    double gradn_1 = (Drag(sz,1)-Drag(sz-2,1))/(Drag(sz,0)-Drag(sz-2,0));

    double curv0 = (grad0_1-grad0)/(Drag(1,0)-Drag(0,0));
    double curvn = (gradn-gradn_1)/(Drag(sz,0)-Drag(sz-1,0));

    if (i==0) {
        // Extrapolation using curvature
        cd = Drag(0,1) + grad0*(alph-Drag(0,0)) + pow(alph-Drag(0,0),2)/2*curv0; 
    } else if (i == Drag.rows()) {
        // Extrapolation using curvature
        cd = Drag(sz,1) + gradn*(alph-Drag(sz,0)) + pow(alph-Drag(sz,0),2)/2*curvn; 
    }else {
        cd = Drag(i-1,1) + (Drag(i,1)-Drag(i-1,1))/(Drag(i,0)-Drag(i-1,0))*(alph-Drag(i-1,0));
    }
    return cd;
}



double inputData::getCDgrad(double a /*Radians*/){

    double grad;
    double step = 0.1/180.0*pi;

    grad = (getCD(a+step)-getCD(a-step))/(2*step);

    return grad;
}

double inputData::getCLgrad(double a /*Radians*/){
    double grad;

    double step = 0.1/180.0*pi;

    grad = (getCL(a+step)-getCL(a-step))/(2*step);

    return grad;
}

double inputData::getAeff(){
    double aeff;
    aeff = -atan((U_inf*sin(alpha+theta)-wv)/(-U_inf*cos(alpha+theta)-uv))-theta;
    return aeff;
}

void inputData::initialize(){
    calcDamp();
    computeEigen();
    calcSteadyLoads();
    // Calculate Partial solution
    partialSol = K.householderQr().solve(Q);
}
