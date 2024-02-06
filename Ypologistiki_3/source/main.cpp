#include "utils.h"
#include "data.h"

using namespace std;

// Define global variables
std::vector<double> x;
std::vector<double> S;
std::vector<double> U1_p;
std::vector<double> U2_p;
std::vector<double> U3_p;
std::vector<double> U1_c;
std::vector<double> U2_c;
std::vector<double> U3_c;
std::vector<double> F1;
std::vector<double> F2;
std::vector<double> F3;
int cells = nodes+1;
std::vector<double> dS;
double dx;


int main()
{
    initialize(); // Initialize variables
    discretize(); // Create grid points
    calcS(); // Calculate values of section area for each node position

    updateBC(); // Assign primitive values to ghost cells
    calcCons(); // Calculate conservative variables for all cells

    roeFlux(); // Calculate face fluxes using Roe's scheme

    std::vector<double> rkConst(4);
    rkConst[0] = 0.1084; rkConst[1] = 0.2602; rkConst[2] = 0.5052; rkConst[3] = 1;

    // Initialize vectors to hold variables of 
    // previous timestep
    std::vector<double> Uold1(cells);
    std::vector<double> Uold2(cells);
    std::vector<double> Uold3(cells);

    for (int it = 1; it <= maxIter; it ++)
    {
        Uold1 = U1_c;
        Uold2 = U2_c;
        Uold3 = U3_c;

        for (int rk = 0; rk<4; rk++)
        {
            // Calculate conservative variables for
            // the next RK step for all (non-ghost) cells
            for (int cl=1; cl<cells-1; cl++)
            {
                U1_c[cl] = Uold1[cl] - dt*rkConst[rk]* (F1[cl]-F1[cl-1])/dx;
                U2_c[cl] = Uold2[cl] + dt*rkConst[rk]* (U3_p[cl]*(dS[cl-1]+dS[cl])/2 -(F2[cl]-F2[cl-1])/dx) ;
                U3_c[cl] = Uold3[cl] - dt*rkConst[rk]* (F3[cl]-F3[cl-1])/dx;
            }
            // Re-calculate values for next RK step
            calcPrim(); 
            updateBC();
            calcCons();
            roeFlux();

        }
        // Handle result writing
        std::cout << it << " " << it*dt << std::endl;
        if (it%writeEvery == 0)
            writeResults(it);
    }

}
