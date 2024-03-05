#include "utilities.h"

void Theodorsen::calcNextY() {
    for (int i = 0; i<2; i++) {
        y[i] = (b[i]*A[i]*2*target->U_inf/target->chord*target->getAeff()+y[i]/dt)/(1/dt+b[i]*2*target->U_inf/target->chord);
    }
}

// Calculate unsteady aerodynamic loads
void Theodorsen::calcLoads() {

    ae = target->getAeff()*(1-A[0]-A[1])+y[0]+y[1];

    double h_dd = target->wa*cos(target->theta) + target->ua*sin(target->theta); // Flapwise acceleration

    double cl = grad0*(ae-a0) - pi*target->chord*h_dd/(2*pow(target->U_inf,2));
    double cd = target->getCD(ae) + cl*(target->getAeff()-ae);

    L = 0.5*target->rho*pow(target->U_inf,2)*cl*target->chord;
    D = 0.5*target->rho*pow(target->U_inf,2)*cd*target->chord;

    target->Q(0) = L*sin(ae+target->theta)-D*cos(ae+target->theta);
    target->Q(1) = L*cos(ae+target->theta)+D*sin(ae+target->theta);

    // File output
    std::ofstream loadsTh("loadsTheodorsen.dat", std::ofstream::app);
    loadsTh << target->t_tot << std::setw(25) << target->Q(0) << std::setw(25) << target->Q(1) << std::endl;
    loadsTh.close();

}

void Theodorsen::calcDamp() {

    // Theodorsen derivation Aeff

    // Aeff partial derivatives
    double partaU = - (target->U_inf*sin(target->alpha+target->theta)-target->wv)/(pow(target->U_inf,2) + pow(target->uv,2) + pow(target->wv,2) + 2*target->U_inf*(cos(target->alpha+target->theta)*target->uv - sin(target->alpha+target->theta)*target->wv));

    double partaW = - (target->U_inf*cos(target->alpha+target->theta)+target->uv)/(pow(target->U_inf,2) + pow(target->uv,2) + pow(target->wv,2) + 2*target->U_inf*(cos(target->alpha+target->theta)*target->uv - sin(target->alpha+target->theta)*target->wv));


    // Cl, Cd partial derivatives
    double partClU = grad0*partaU;
    double partClW = grad0*partaW;

    double partCdU = target->getCDgrad(ae)*partaU + partClU*(target->getAeff()-ae);
    double partCdW = target->getCDgrad(ae)*partaW + partClW*(target->getAeff()-ae);

    // Lift, Drag derivatives
    double partLU = 0.5*target->rho*pow(target->U_inf,2)*target->chord*partClU;
    double partLW = 0.5*target->rho*pow(target->U_inf,2)*target->chord*partClW;

    double partDU = 0.5*target->rho*pow(target->U_inf,2)*target->chord*partCdU;
    double partDW = 0.5*target->rho*pow(target->U_inf,2)*target->chord*partCdW;

    // C matrix assembly
    target->C(0,0) = partLU*sin(ae+target->theta) + L*cos(ae+target->theta)*partaU 
                    - partDU*cos(ae+target->theta) + D*sin(ae+target->theta)*partaU;

    target->C(0,1) = partLW*sin(ae+target->theta) + L*cos(ae+target->theta)*partaW 
                    - partDW*cos(ae+target->theta) + D*sin(ae+target->theta)*partaW;

    target->C(1,0) = partLU*cos(ae+target->theta) - L*sin(ae+target->theta)*partaU 
                    + partDU*sin(ae+target->theta) + D*cos(ae+target->theta)*partaU;

    target->C(1,1) = partLW*cos(ae+target->theta) - L*sin(ae+target->theta)*partaW 
                    + partDW*sin(ae+target->theta) + D*cos(ae+target->theta)*partaW;

    target->C = -target->C;

    // File output
    std::ofstream dampTh("dampingTheodorsen.dat", std::ofstream::app); 
    dampTh << target->t_tot << std::setw(25) << target->C(0,0) << std::setw(25) << target->C(0,1) << std::setw(25) << target->C(1,0) << std::setw(25) << target->C(1,1) << std::endl;
    dampTh.close();

}
