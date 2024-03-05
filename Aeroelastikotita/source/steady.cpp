#include "utilities.h"

// Calculates steady damping matrix
void inputData::calcDamp() {
    double partaU, partaW, aeff;
    
    // Aeff partial derivatives

    partaU = - (U_inf*sin(alpha+theta)-wv)/(pow(U_inf,2) + pow(uv,2) + pow(wv,2) + 2*U_inf*(cos(alpha+theta)*uv - sin(alpha+theta)*wv));
    partaW = - (U_inf*cos(alpha+theta)+uv)/(pow(U_inf,2) + pow(uv,2) + pow(wv,2) + 2*U_inf*(cos(alpha+theta)*uv - sin(alpha+theta)*wv));
    aeff = -atan((U_inf*sin(alpha+theta)-wv)/(-U_inf*cos(alpha+theta)-uv))-theta;

    // Construction of C Matrix
    C(0,0) = -0.5*rho*chord*pow(U_inf,2)*partaU * (sin(aeff+theta)*(getCLgrad(aeff)+getCD(aeff)) + cos(aeff+theta)*(getCL(aeff)-getCDgrad(aeff)) );

    C(0,1) = -0.5*rho*chord*pow(U_inf,2)*partaW * (sin(aeff+theta)*(getCLgrad(aeff)+getCD(aeff)) + cos(aeff+theta)*(getCL(aeff)-getCDgrad(aeff)) );

    C(1,0) = -0.5*rho*chord*pow(U_inf,2)*partaU * (sin(aeff+theta)*(getCDgrad(aeff)-getCL(aeff)) + cos(aeff+theta)*(getCD(aeff)+getCLgrad(aeff)) );

    C(1,1) = -0.5*rho*chord*pow(U_inf,2)*partaW * (sin(aeff+theta)*(getCDgrad(aeff)-getCL(aeff)) + cos(aeff+theta)*(getCD(aeff)+getCLgrad(aeff)) );

    // File output
    std::ofstream dampSt("dampingSteady.dat", std::ofstream::app); 
    dampSt << t_tot << std::setw(25) << C(0,0) << std::setw(25) << C(0,1) << std::setw(25) << C(1,0) << std::setw(25) << C(1,1) << std::endl;
    dampSt.close();
}

// Calculates steady aerodynamic loads
void inputData::calcSteadyLoads() {

    double aeff = getAeff();
    double L, D;
    L = 0.5*getCL(aeff)*rho*pow(U_inf,2)*chord;
    D = 0.5*getCD(aeff)*rho*pow(U_inf,2)*chord;

    Q(0) = L*sin(aeff+theta)-D*cos(aeff+theta);
    Q(1) = L*cos(aeff+theta)+D*sin(aeff+theta);

    // File output

    std::ofstream loadsSt("loadsSteady.dat", std::ofstream::app);
    loadsSt << t_tot << std::setw(25) << Q(0) << std::setw(25) << Q(1) << std::endl;
    loadsSt.close();

}
