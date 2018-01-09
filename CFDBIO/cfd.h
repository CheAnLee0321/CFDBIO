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
    void CFD_NSInitialGuessComplex3D();
    void CFD_NewAndInitialize();
    void CFD_SIMPLE2D();
    void CFD_SIMPLER2D();
    void CFD_SIMPLER3D();

    //Vx
    double CFD_VxSolver2D();
    double CFD_VxGaussSeidel2D();
    double CFD_VxGaussSeidelInner2D(int i, int j);
    double CFD_VxSolver3D();
    double CFD_VxGaussSeidel3D();
    double CFD_VxGaussSeidelInner3D(int i, int j, int k);
    void   CFD_dVxCalculation2D();
    void   CFD_dVxCalculation3D();
    void   CFD_pseudoVx2D();
    void   CFD_pseudoVx3D();

    //Vy
    double CFD_VySolver2D();
    double CFD_VyGaussSeidel2D();
    double CFD_VyGaussSeidelInner2D(int i, int j);
    double CFD_VySolver3D();
    double CFD_VyGaussSeidel3D();
    double CFD_VyGaussSeidelInner3D(int i, int j, int k);
    void   CFD_dVyCalculation2D();
    void   CFD_dVyCalculation3D();
    void   CFD_pseudoVy2D();
    void   CFD_pseudoVy3D();

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
    double CFD_dPSolver3D();
    double CFD_dPGaussSeidel3D();
    double CFD_dPGaussSeidelInner3D(int i, int j, int k);

    //P
    double CFD_PSolver2D();
    double CFD_PGaussSeidel2D();
    double CFD_PGaussSeidelInner2D(int i, int j);
    double CFD_PSolver3D();
    double CFD_PGaussSeidel3D();
    double CFD_PGaussSeidelInner3D(int i, int j, int k);

    //Correction
    void CFD_VxCorrection2D();
    void CFD_VyCorrection2D();
    void CFD_PCorrection2D();
    void CFD_VxCorrection3D();
    void CFD_VyCorrection3D();
    void CFD_VzCorrection3D();

    //Calculation
    void CFD_VelocityCalculation2D();
    void CFD_VelocityCalculation3D();

    //I/O
    void CFD_PrintVx2D(const char *path);
    void CFD_PrintVy2D(const char *path);
    void CFD_PrintMaterialComplex2D(const char *path);
    void CFD_PrintMaterialComplex3D(const char *path);
    void CFD_PrintVelocity2D(const char *path);
    void CFD_PrintVelocity3D(const char *path);
    void CFD_PrintPotential2D(const char *path);
    void CFD_PrintPotential3D(const char *path);
    void CFD_ReadMaterialComplex2D(const char *path);
    void CFD_ReadMaterialComplex3D(const char *path);
    void CFD_ReadVelocity3D(const char *path);
    void CFD_ReadPotential3D(const char *path);
    void CFD_ReadVelocity2D(const char *path);
    void CFD_ReadPotential2D(const char *path);
    void CFD_MaxSpeed2D();
    void CFD_MaxSpeed3D();

    //CFD Poisson Solver
    void CFD_Initialize();
    void CFD_ACPoissonInitialGuess2D();
    void CFD_ACPoissonInitialGuess3D();
    double CFD_PoissonSolverComplex2D();
    double CFD_PoissonSolverComplex3D();
    double CFD_PoissonGaussSeidelComplex2D();
    double CFD_PoissonGaussSeidelComplex3D();
    complex<double> CFD_PoissonGaussSeidelInnerComplex2D(int i, int j);
    complex<double> CFD_PoissonGaussSeidelInnerComplex3D(int i, int j, int k);

    void ACPoissonBC2D();
    void ACPoissonBC3D();

    int ACDCPhaseFlag; // DC=1, AC=2

private:

    //Boundary
    void CFD_VxBoundary2D();
    void CFD_VyBoundary2D();
    void CFD_VxBoundaryDC2D();
    void CFD_VyBoundaryDC2D();
    void CFD_VxBoundaryAC2D();
    void CFD_VyBoundaryAC2D();
    void CFD_VxBoundary3D();
    void CFD_VyBoundary3D();
    void CFD_VzBoundary3D();

    //coef
    double CFD_uAcoef2D(int i, int j);
    double CFD_vAcoef2D(int i, int j);
    double CFD_BPcoef2D(int i, int j);
    double CFD_Bcoef2D(int i, int j);
    double CFD_uAcoef3D(int i, int j, int k);
    double CFD_vAcoef3D(int i, int j, int k);
    double CFD_wAcoef3D(int i, int j, int k);
    double CFD_BPcoef3D(int i, int j, int k);
    double CFD_Bcoef3D(int i, int j, int k);

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
    double Vac, CapLambda;
    double Lambda_D1, Lambda_D2, D_KCl, Omega, CDL, COxide, CTotal, Conductivity,CC0;
};

#endif // CFD_H
