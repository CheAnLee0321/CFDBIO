#ifndef CFD_H
#define CFD_H
#include <complex>
#include "fvmesh.h"

struct Electrolyte {
    double k;   // permitivitty
    complex <double> phi; // potential, Ei
    double u;   // quasi-fermi level, exp(-Ef/Vt)
    double v;   //                  , exp(Ef/Vt)
    double r;   // recombination rate, SRH
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

    /*
    // Basic funstions
    void NSParameter();
    void NSInitialGuessComplex2D();
    void NSInitialGuessComplex2D3E();
    void NSInitialGuess3D();
    void NSInitialGuessComplex3D_CompoundEyes();
    void NSInitialGuessComplex3D_Gemini();
    void NSInitialGuessComplex3D_SuperGuardRing();
    void SIMPLE();
    void SIMPLE2D();
    void SIMPLE3D();
    void SIMPLER();
    void SIMPLER2D();
    void SIMPLER3D();

    //Vx
    double VxSolver2D();
    double VxGaussSeidel2D();
    double VxGaussSeidelInner2D(int i, int j);
    void   dVxCalculation2D();
    void   pseudoVx2D();

    double VxSolver3D();
    double VxGaussSeidel3D();
    double VxGaussSeidelInner3D(int i, int j, int k);
    void   dVxCalculation3D();
    void   pseudoVx3D();

    //Process
    //Vy
    double VySolver2D();
    double VyGaussSeidel2D();
    double VyGaussSeidelInner2D(int i, int j);
    void   dVyCalculation2D();
    void   pseudoVy2D();

    double VySolver3D();
    double VyGaussSeidel3D();
    double VyGaussSeidelInner3D(int i, int j, int k);
    void   dVyCalculation3D();
    void   pseudoVy3D();

    //Vz
    double VzSolver3D();
    double VzGaussSeidel3D();
    double VzGaussSeidelInner3D(int i, int j, int k);
    void   dVzCalculation3D();
    void   pseudoVz3D();

    //dP
    double dPSolver2D();
    double dPGaussSeidel2D();
    double dPGaussSeidelInner2D(int i, int j);


    double dPSolver3D();
    double dPGaussSeidel3D();
    double dPGaussSeidelInner3D(int i, int j, int k);

    //P
    double PSolver2D();
    double PGaussSeidel2D();
    double PGaussSeidelInner2D(int i, int j);


    double PSolver3D();
    double PGaussSeidel3D();
    double PGaussSeidelInner3D(int i, int j, int k);

    //coef
    double uAcoef2D(int i, int j);
    double vAcoef2D(int i, int j);
    double BPcoef2D(int i, int j);
    double Bcoef2D(int i, int j);

    double uAcoef3D(int i, int j, int k);
    double vAcoef3D(int i, int j, int k);
    double wAcoef3D(int i, int j, int k);
    double BPcoef3D(int i, int j, int k);
    double Bcoef3D(int i, int j, int k);

    //Correction
    void VxCorrection2D();
    void VyCorrection2D();
    void PCorrection2D();
    void VxCorrection3D();
    void VyCorrection3D();
    void VzCorrection3D();
    void PCorrection3D();
    void Correction3D();

    //Boundary
    void VxBoundary2D();
    void VyBoundary2D();
    void VxBoundaryDC2D();
    void VyBoundaryDC2D();
    void VxBoundaryAC2D();
    void VyBoundaryAC2D();
    void VxBoundary3D();
    void VyBoundary3D();
    void VzBoundary3D();
    void VxBoundaryDC3D();
    void VyBoundaryDC3D();
    void VzBoundaryDC3D();
    void VxBoundaryAC3D();
    void VyBoundaryAC3D();
    void VzBoundaryAC3D();

protected:


    int NSloop, MaxIter;
    double SimTolVelocity, SimTolPressure, Volt, Zeta1, Zeta2;
    double alphadP, alphaP, alphaVx, alphaVy, alphaVz;
    double alphaPc, alphaPhi;
    double Pbreak,dPbreak,Vxbreak,Vybreak,Vzbreak;

    */

protected:

    Electrolyte *solution;

    double SimTolVelocity, SimTolPressure;
    double alphadP, alphaP, alphaVx, alphaVy, alphaVz,alphaPc, alphaPhi;
    double Pbreak,dPbreak,Vxbreak,Vybreak,Vzbreak;

};

#endif // CFD_H
