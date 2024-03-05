#include "utilities.h"


void inputData::computeEigen(){

    // Initialize help matrices
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> R;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> newR;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> invL;

    R.resize(2*siz, 2*siz);
    newR.resize(2*siz, 2*siz);
    L.resize(2*siz, 2*siz);
    invL.resize(2*siz, 2*siz);

    if(theodor)
    {
        C.setZero();
    }

    // Construct A Matrix
    R.setZero(); L.setIdentity();
    R.block(0,0,2,2) = -C;
    R.block(0,2,2,2) = -K;
    R(2,0) = 1;
    R(3,1) = 1;

    L(0,0) = M(0,0); L(1,1) = M(1,1); 
    invL = L.colPivHouseholderQr().inverse();
    
    newR = invL*R;

    // Compute eigenvalues of A
    Eigen::EigenSolver<Eigen::MatrixXd> ces;
    ces.setMaxIterations(500);
    ces.compute(newR, true);

    // Save eigenvectors in data class
    for (int i =0;i<4;i++)
    {
        eigenvectors(i,0).set(ces.eigenvectors()(i,0));
        eigenvectors(i,1).set(ces.eigenvectors()(i,2));
    }

    eigenvalues[0] = ces.eigenvalues()[0];
    eigenvalues[1] = ces.eigenvalues()[2];

}

// Solves for the homogenous solution coefficients
// with boundary conditions as input
void inputData::computeConstants() {

    Eigen::Matrix<double,4,4> coeff;
    Eigen::Matrix<double,4,1> rhs;
    Eigen::Matrix<double,4,1> temp;

    coeff.setZero();

    // Solves x=x_0 for t0 system of equations
    for (int i = 0;i<4;i++)
    {
        coeff(i,0) = eigenvectors(i,0).len*cos(eigenvectors(i,0).theta);
        coeff(i,1) = eigenvectors(i,0).len*sin(eigenvectors(i,0).theta);
        coeff(i,2) = eigenvectors(i,1).len*cos(eigenvectors(i,1).theta);
        coeff(i,3) = eigenvectors(i,1).len*sin(eigenvectors(i,1).theta);
    }

    // Set boundary conditions as rhs
    rhs(0,0) = uv;
    rhs(1,0) = wv;
    rhs(2,0) = u;
    rhs(3,0) = w;
    
    // Calculate Partial solution
    if (theodor)
    {
        partialSol = K.householderQr().solve(Q);
    }
    
    // Add partial solution to the rhs
    rhs(2,0) -= partialSol(0);
    rhs(3,0) -= partialSol(1);

    // Solve system -> obtain coefficients
    temp = coeff.householderQr().solve(rhs);
    A1 = temp(0,0);
    B1 = temp(1,0);
    A2 = temp(2,0);
    B2 = temp(3,0);

    // Output state data for plotting
    std::ofstream out("state.dat");

    double ang, x1, x2, y1, y2;

    for (double t = 0; t<360; t ++)
    {
        ang = t*pi/180;
        x1 = eigenvectors(2,0).len*(A1*cos(ang+eigenvectors(2,0).theta) + B1*sin(ang+eigenvectors(2,0).theta));
        x2 = eigenvectors(2,1).len*(A2*cos(ang+eigenvectors(2,1).theta) + B2*sin(ang+eigenvectors(2,1).theta));
        y1 = eigenvectors(3,0).len*(A1*cos(ang+eigenvectors(3,0).theta) + B1*sin(ang+eigenvectors(3,0).theta));
        y2 = eigenvectors(3,1).len*(A2*cos(ang+eigenvectors(3,1).theta) + B2*sin(ang+eigenvectors(3,1).theta));
        out << x1 << " " << y1 << " " <<  x2 << " " << y2 << std::endl;
    }
    out.close();
}
// Newmark integration (Unused)
void inputData::calcNextStep(){
    
    // Newmark integration method
    
    // Newmark constants
    double beta = 1.0/4.0;
    double Gamma = 0.5;

    Eigen::Matrix<double, 2, 2> Keff;
    Eigen::Vector<double, 2> Qeff, x, xv, temp;

    t_tot += dt;
 
    // Intermediate predictor values
    x(0) = u + dt*uv + (0.5-beta)*pow(dt,2)*ua;
    x(1) = w + dt*wv + (0.5-beta)*pow(dt,2)*wa;

    xv(0) = uv + (1-Gamma)*dt*ua;
    xv(1) = wv + (1-Gamma)*dt*wa;

    // Construct Effective stiffness and load matrices
    Keff = M/(beta*pow(dt,2)) + C*Gamma/(beta*dt) + K;

    Qeff = Q + (M/(beta*pow(dt,2)) + Gamma*C/(beta*dt))*x - C*xv;

    temp = Keff.householderQr().solve(Qeff);

    u = temp(0);
    w = temp(1);

    // Update velocity and accel values
    uv = xv(0) + Gamma/(beta*dt)*(u-x(0));
    wv = xv(1) + Gamma/(beta*dt)*(w-x(1));

    ua = (u-x(0))/(beta*dt);
    wa = (w-x(1))/(beta*dt);


    std::ofstream out("time.dat", std::ofstream::app); 

     out << t_tot << std::setw(25) << u << std::setw(25) << w << std::setw(25) << uv << std::setw(25) << wv << std::setw(25) << ua << std::setw(25) << wa << std::endl;
    out.close();

}

void inputData::nextAnalyticalStep(std::string filename) {

    std::vector<double> var(4);

    t_tot += dt;

    // Calculate homogenous solution for each variable
    for (int i = 0; i<4; i++)
    {
        var[i] = eigenvectors(i,0).len*std::exp(eigenvalues[0].real()*dt)*(A1*cos(eigenvalues[0].imag()*dt+eigenvectors(i,0).theta)
                + B1*sin(eigenvalues[0].imag()*dt + eigenvectors(i,0).theta))

               + eigenvectors(i,1).len*std::exp(eigenvalues[1].real()*dt)*(A2*cos(eigenvalues[1].imag()*dt+eigenvectors(i,1).theta)
               + B2*sin(eigenvalues[1].imag()*dt + eigenvectors(i,1).theta));
                
    }


    // Add partial solution
    var[2] += partialSol(0);
    var[3] += partialSol(1);

    // Calculate acceleration values (FD)
    ua = (var[0]-uv)/dt;
    wa = (var[1]-wv)/dt;

    // Update solution state variables
    u = var[2]; w = var[3]; uv = var[0]; wv = var[1]; 

    // File output
    std::ofstream outAn(filename, std::ofstream::app); 

    outAn << t_tot << std::setw(25) << u << std::setw(25) << w << " " << uv << " " << wv << " " << ua << " " << wa << std::endl;

    outAn.close();
}
