#include "utils.h"
#include "data.h"
#include <fstream>
#include <string>

// Declare global variables
extern std::vector<double> x;
extern std::vector<double> S;
extern std::vector<double> U1_p;
extern std::vector<double> U2_p;
extern std::vector<double> U3_p;
extern std::vector<double> U1_c;
extern std::vector<double> U2_c;
extern std::vector<double> U3_c;
extern std::vector<double> F1;
extern std::vector<double> F2;
extern std::vector<double> F3;
extern int cells;
extern std::vector<double> dS;
extern double dx;

void initialize()
{
    // Allocate memory for all vectors
    x.resize(nodes);
    S.resize(nodes);
    dS.resize(nodes);

    F1.resize(nodes);
    F2.resize(nodes);
    F3.resize(nodes);

    U1_p.resize(cells);
    U2_p.resize(cells);
    U3_p.resize(cells);
    U1_c.resize(cells);
    U2_c.resize(cells);
    U3_c.resize(cells);

    // Initialize all vectors with their respective values
    std::fill(U1_p.begin(), U1_p.end(), densOut);
    std::fill(U2_p.begin(), U2_p.end(), 0);
    std::fill(U3_p.begin(), U3_p.end(), presOut);

}
void discretize() {
    dx = tubeLength/(nodes-1);
    for (int i =0; i<nodes; i++)
    {
        x[i] = i*dx;
    }
}

void calcS()
{
    double d = 0;
    double x_s = 0;
    int size = (int) S.size();
    std::ofstream sOut("results/S.dat");

    for (int i = 0; i<size; i ++)
    {
        x_s = x[i] - cc;
        // Calculate section area
        S[i] = k + a*pow(x_s,2)*(pow(x_s,2)-pow(b,2));
        // Calculate diameter
        d = sqrt(4*S[i]/pi);
        // Save area derivative
        dS[i] = 2 * a * x_s *( pow(x_s,2)-pow(b,2) ) + 2 * a * pow(x_s,3);
        // File output
        sOut << x[i] << " " << d << " " << S[i] << std::endl;
    }
    sOut.close();
}

void calcCons()
{
    
    double Sm = 1; // Mean area of cell
    
    for (int i = 0; i < cells; i ++)
    {
        // Interpolate area of adjacent nodes 
        // as the mean area of the cell
        if (i!=0 && i!= cells-1)
        {
            Sm = 0.5 * (S[i-1]+S[i]);
        }else if (i==0) // If first ghost cell take the area of the first node
        {
            Sm = S[0];
        }else // If last ghost -> area of last node
        {
            Sm = S[nodes];
        }

        U1_c[i] = U1_p[i]*Sm;
        U2_c[i] = U2_p[i]*U1_p[i]*Sm;
        U3_c[i] = ( U3_p[i]/Gamma1 + 0.5 * U1_p[i] * pow( U2_p[i], 2) ) * Sm;
    }

}

void calcPrim()
{
    double Sm = 1; // Mean area of cell
    
    for (int i = 0; i < cells; i ++)
    {
        // Interpolate area of adjacent nodes 
        // as the mean area of the cell
        if (i!=0 && i!= cells-1)
        {
            Sm = 0.5 * (S[i-1]+S[i]);
        }else if (i==0)
        {
            Sm = S[0];
        }else
        {
            Sm = S[nodes];
        }

        U1_p[i] = U1_c[i]/Sm;
        U2_p[i] = U2_c[i]/U1_c[i];
        U3_p[i] = Gamma1 * ( U3_c[i] - 0.5*pow(U2_c[i] ,2) / U1_c[i] )/Sm;
    }

}

void calcFlux()
{
    std::vector<double> FL(3), FR(3);
    double VN_L, VN_R, NORM_L, NORM_R;
    NORM_L = 1;
    NORM_R = 1;
    for (int i=0; i<nodes; i++)
    {
        VN_L = NORM_L*U2_p[i];
        VN_R = NORM_R*U2_p[i+1];
        FL[0] = U1_p[i]*VN_L*S[i];
        FL[1] = (U1_p[i] * U2_p[i] * VN_L + NORM_L*U3_p[i])*S[i];
        FL[2] = (U3_c[i]/S[i] + U3_p[i]) * VN_L * S[i];

        FR[0] = U1_p[i+1]*VN_R*S[i];
        FR[1] = (U1_p[i+1] * U2_p[i+1] * VN_R + NORM_R*U3_p[i+1])*S[i];
        FR[2] = (U3_c[i+1]/S[i] + U3_p[i+1]) * VN_R * S[i];

        F1[i] = 0.5*(FL[0]+FR[0]);
        F2[i] = 0.5*(FL[1]+FR[1]);
        F3[i] = 0.5*(FL[2]+FR[2]);
    }
}

// Assemble left matrix
Eigen::Matrix<double, 3, 3> asmLeft(double U, double dens, double C)
{
    Eigen::Matrix<double, 3, 3> left;
    
    left(0,0) = 1 - 0.5 * Gamma1 * pow(U,2)/pow(C,2);
    left(0,1) = Gamma1 * U / pow(C,2);
    left(0,2) = - Gamma1 / pow(C,2);

    left(1,0) = (0.5 * Gamma1 * pow(U,2) - U*C) * 1/(dens*C);
    left(1,1) = (C - Gamma1 * U)/(dens*C);
    left(1,2) = Gamma1/(dens*C);

    left(2,0) = - (0.5 * Gamma1 * pow(U,2) + U*C)/(dens*C);
    left(2,1) = (C+Gamma1*U)/(dens*C);
    left(2,2) = - Gamma1/(dens*C);

    return left;
}

// Assemble right matrix
Eigen::Matrix<double, 3, 3> asmRight(double U, double dens, double C)
{

    Eigen::Matrix<double, 3, 3> right;

    right(0,0) = 1;
    right(0,1) = 0.5*dens/C;
    right(0,2) = -0.5*dens/C;

    right(1,0) = U;
    right(1,1) = 0.5*(U+C)*dens/C;
    right(1,2) = -0.5*(U-C)*dens/C;

    
    right(2,0) = 0.5*pow(U,2);
    right(2,1) = (0.5 * pow(U,2) + U*C + pow(C,2)/Gamma1)*0.5*dens/C;
    right(2,2) = -(0.5 * pow(U,2) - U*C + pow(C,2)/Gamma1)*0.5*dens/C;

    return right;
}

void roeFlux()
{
    calcFlux();
    double r,u,HL,HR,H,c,delta;
    Eigen::Matrix<double, 3, 3> right, left, lambda, alpha;
    Eigen::Matrix<double, 3, 1> ULR;
    Eigen::Matrix<double, 3, 1> roeCor;
    right.setZero(); left.setZero(); lambda.setZero(); ULR.setZero(); alpha.setZero();

    for (int i = 0; i<nodes; i++)
    {
        r = pow(U1_p[i],0.5)*pow(U1_p[i+1],0.5);
        u = ( pow(U1_p[i],0.5)*U2_p[i] + pow(U1_p[i+1],0.5)*U2_p[i+1] )
            /(pow(U1_p[i],0.5) + pow(U1_p[i+1],0.5));
        HL = (U3_c[i]/S[i]+U3_p[i])/U1_p[i]; 
        HR = (U3_c[i+1]/S[i]+U3_p[i+1])/U1_p[i+1]; 
        H = (pow(U1_p[i],0.5)*HL + pow(U1_p[i+1],0.5)*HR )
            /(pow(U1_p[i],0.5) + pow(U1_p[i+1],0.5));
        c = sqrt(Gamma1*(H-0.5*pow(u,2)));
        

        delta = 0.05*c;
        
        // Assemble right and left matrices
        left = asmLeft(u, r, c);
        right = asmRight(u, r, c);
        
        // Perform corrections to handle shock waves
        lambda(0,0) = abs(u);
        if (lambda(0,0)<=delta)
            lambda(0,0) = (pow(lambda(0,0),2)+pow(delta,2))/(2*delta);

        lambda(1,1) = abs(u+c);
        if (lambda(1,1)<=delta)
            lambda(1,1) = (pow(lambda(1,1),2)+pow(delta,2))/(2*delta);

        lambda(2,2) = abs(u-c);
        if (lambda(2,2)<=delta)
            lambda(2,2) = (pow(lambda(2,2),2)+pow(delta,2))/(2*delta);

        // Difference of right - left conservative vars vector
        ULR(0,0) = U1_c[i+1]-U1_c[i];
        ULR(1,0) = U2_c[i+1]-U2_c[i];
        ULR(2,0) = U3_c[i+1]-U3_c[i];

        // Calc A matrix
        alpha = right*lambda*left;
        roeCor = alpha*ULR;

        // Calculate final fluxes for face
        F1[i] = F1[i] - 0.5*roeCor(0,0);
        F2[i] = F2[i] - 0.5*roeCor(1,0);
        F3[i] = F3[i] - 0.5*roeCor(2,0);
        
    }

}

void updateBC()
{
    // Inflow BC
    U1_p[0] = densIn; 
    U2_p[0] = U2_p[1]; 
    U3_p[0] = presIn; 
    // Outflow BC
    U1_p.back() = U1_p[cells-2];
    U2_p.back() = U2_p[cells-2];
    U3_p.back() = presOut;
}


void writeResults(int iter)
{
    std::string primName = "primitiveVars";
    std::string extension = ".dat";
    std::string filepath = "results/";
    std::string prFileName = filepath.append(primName).append(std::to_string(iter)).append(extension);
    std::ofstream outVar(prFileName);


    outVar << "Time Iteration: " << iter <<std::endl;
    outVar << "Time: " << (iter)*dt <<std::endl;
    outVar << std::setw(25) << "Position"
            << std::setw(30) << "Density"
            << std::setw(30) << "Velocity"
            << std::setw(30) << "Speed of sound"
            << std::setw(30) << "Mach"
            << std::setw(30) << "Pressure" << std::endl;

    double xm;
    for (int j=0; j<cells; j++)
    {
        if (j!=0 && j!=cells-1)
        {
            xm = (x[j-1]+x[j])/2;
        }else if (j==0)
        {
            xm = -(x[0]+x[1])/2;
        }else
            xm = -(x[nodes-1]+x[nodes-1])/2 + x[nodes-1];

        double u_c = sqrt(Gamma*U3_p[j]/(U1_p[j]));
        outVar << std::setw(25) << xm
            << std::setw(30) << U1_p[j]
            << std::setw(30) << U2_p[j]
            << std::setw(30) << u_c // Sound speed
            << std::setw(30) << U2_p[j]/u_c // Mach number
            << std::setw(30) << U3_p[j] << std::endl;
    }
    outVar.close();
}


