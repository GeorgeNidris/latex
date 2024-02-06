#pragma once

// Physical properties
#define Gamma 1.4
#define Gamma1 0.4
#define R 287 // J/(kgK)
#define presIn 400000.0 //4 bar to Pa
#define presOut 330000.0 // Pa
#define tempIn 286.0 // K
#define densIn presIn/(R*tempIn) // kg/m3
#define densOut presOut/(R*tempIn) // kg/m3

// Simulation parameters
#define nodes 801
#define dt 0.000001
#define writeEvery 100
#define maxIter 45000

// Geometry parameters
#define tubeLength 2.0
#define a 0.216
#define b 1.03
#define cc 0.97
#define k 0.09
