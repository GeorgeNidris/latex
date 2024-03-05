#include "utilities.h"
#include <complex>
#include <fstream>
#include <ostream>


typedef std::complex<double> cmp;

int main(int argc, char *argv[]) {

    double time = std::stod(argv[1]);
    int sz;

    sz = 2;
    // Calculation classes construction
    inputData data(sz);
    inputData alt(sz);

    std::string filename = "matrix.dat";


    // Read inputs
    readMat(sz, data.M, data.K, filename);
    readMat(sz, alt.M, alt.K, filename);

    data.readInput("input.dat");
    alt.readInput("input.dat");

    // Read aerodynamic data
    std::string clpath = "doc/CL.dat";
    std::string cdpath = "doc/CD.dat";
    readAero(clpath, cdpath, data);
    readAero(clpath, cdpath, alt);

    // Initialize solutions
    // Calculates Damping, loads
    // and steady state solution
    data.initialize();
    alt.initialize();
    data.theodor = true; // data class will use Theodorsen equations

    // Δημιουργία και αρχικοποίηση αρχείων εξόδου

        std::ofstream outE("theodorsen.dat"); 
        outE << "Time (s) " << std::setw(25) << "U" << std::setw(25) << "W" << std::endl;
        outE.close();

        std::ofstream outA("eigenSteady.dat"); 
        outA << "Time (s) " << std::setw(25) << "U" << std::setw(25) << "W" << std::endl;
        outA.close();
    
    // Unsteady class calculation initialization
    Theodorsen theodorsen(&data);
    
    // Main time integration loop
    while (data.t_tot < time) 
    {
        std::cout << data.t_tot << "s" << std::setw(25) << data.eigenvalues[0] << std::setw(25) << data.eigenvalues[1] << std::endl;

        theodorsen.calcNextY(); // Calculate next discrete y1, y2
        theodorsen.calcLoads(); 
        data.computeEigen();
        data.computeConstants(); // Compute homogenous solution coefficients
        data.nextAnalyticalStep("theodorsen.dat");

        alt.computeConstants();
        alt.nextAnalyticalStep("eigenSteady.dat");
    }
    
    
    // Angle of attack vs Damping

    // std::ofstream al("damping.dat");
    // std::ofstream freq("frequency.dat");
    // eigen.alpha = -pi/4;
    // while (eigen.alpha < pi/4) {
    //     eigen.calcDamp();
    //     eigen.computeEigen();
    //
    //     al << eigen.alpha/pi*180 << " " << eigen.eigenvalues[0].real() << " " << eigen.eigenvalues[1].real() << std::endl;
    //     freq << eigen.alpha/pi*180 << " " << eigen.eigenvalues[0].imag() << " " << eigen.eigenvalues[1].imag() << std::endl;
    //
    //     eigen.alpha += 2*pi/180;
    //
    // } 
    // al.close();
    // freq.close();

}
