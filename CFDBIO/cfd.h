#ifndef CFD_H
#define CFD_H
#include <complex>
#include "fvmesh.h"

struct CFDElectrolyte {
    double k;   // permitivitty
    complex <double> phi; // potential, Ei
    double u;   // quasi-fermi level, exp(-Ef/Vt)
    double v;   //                  , exp(Ef/Vt)
    double rho; // charge density, C/cm3
    double Crho; // charge density, C/cm3
    double mun; // mobility, nm^2/v.s
    double mup;
    double Ex;  // electric field
    double Ey;
    double Ez;
    double p;   //pressure
    double dp;
    double ui;  //velocity
    double vi;
    double wi;
    double flag;// 1=channel, 2=insulator, 3=electrolyte, 4=sub, 5=analyte;
};

struct Vx {
    double u;
    double du;
    double cu;
    double pu;
};

struct Vy {
    double v;
    double dv;
    double cv;
    double pv;
};

struct Vz {
    double w;
    double dw;
    double cw;
    double pw;
};

class CFD : public FVMesh{

public:

    CFD();
    ~CFD();

    // Basic funstions
    void CFD_NSParameter();
    void CFD_NSInitialGuessComplex2D();
    void CFD_SIMPLE2D();
    void CFD_SIMPLER2D();

    //Vx
    double CFD_VxSolver2D();
    double CFD_VxGaussSeidel2D();
    double CFD_VxGaussSeidelInner2D(int i, int j);
    void   CFD_dVxCalculation2D();
    void   CFD_pseudoVx2D();

    //Vy
    double CFD_VySolver2D();
    double CFD_VyGaussSeidel2D();
    double CFD_VyGaussSeidelInner2D(int i, int j);
    void   CFD_dVyCalculation2D();
    void   CFD_pseudoVy2D();

    //Vz
    double CFD_VzSolver3D();
    double CFD_VzGaussSeidel3D();
    double CFD_VzGaussSeidelInner3D(int i, int j, int k);
    void   CFD_dVzCalculation3D();
    void   CFD_pseudoVz3D();

    //dP
    double CFD_dPSolver2D();
    double CFD_dPGaussSeidel2D();
    double CFD_dPGaussSeidelInner2D(int i, int j);

    //P
    double CFD_PSolver2D();
    double CFD_PGaussSeidel2D();
    double CFD_PGaussSeidelInner2D(int i, int j);

    //coef
    double CFD_uAcoef2D(int i, int j);
    double CFD_vAcoef2D(int i, int j);
    double CFD_BPcoef2D(int i, int j);
    double CFD_Bcoef2D(int i, int j);

    //Correction
    void CFD_VxCorrection2D();
    void CFD_VyCorrection2D();
    void CFD_PCorrection2D();

    //Boundary
    void CFD_VxBoundary2D();
    void CFD_VyBoundary2D();
    void CFD_VxBoundaryDC2D();
    void CFD_VyBoundaryDC2D();
    void CFD_VxBoundaryAC2D();
    void CFD_VyBoundaryAC2D();

    //Calculation
    void CFD_VelocityCalculation2D();

    //I/O
    void CFD_PrintVx2D(const char *path);
    void CFD_PrintVy2D(const char *path);
    void CFD_PrintMaterialComplex2D(const char *path);

    //CFD Poisson Solver
    void CFD_ACPoissonInitialGuess2D();
    double CFD_PoissonSolverComplex2D();
    double CFD_PoissonGaussSeidelComplex2D();
    complex<double> CFD_PoissonGaussSeidelInnerComplex2D(int i, int j);

    void ACPoissonBC2D();

    int ACDCPhaseFlag; // DC=1, AC=2

protected:

    CFDElectrolyte *CFDmaterial;
    Vx *u;
    Vy *v;
    Vz *w;

    int NSloop, MaxIter;
    double SimTolVelocity, SimTolPressure,SimTolPoisson;
    double alphadP, alphaP, alphaVx, alphaVy, alphaVz,alphaPc, alphaPhi;
    double Pbreak,dPbreak,Vxbreak,Vybreak,Vzbreak;

    //ACEO
    double ACFreq, CharFreq;
    double Vac, ElectrodeGap, ElectrodeWidth, CapLambda;
    double Lambda_D1, Lambda_D2, D_KCl, Omega, CDL, COxide, CTotal, Conductivity,CC0;
};

#endif // CFD_H
