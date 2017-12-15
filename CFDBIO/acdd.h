#ifndef ACDD_H
#define ACDD_H
#include "fvmesh.h"
#include <complex>

struct Semiconductor_complex {
    double k;   // permitivitty
    double dop; // doping, -Nai Ndi
    complex <double> phi; // potential, Ei
    complex <double> u;   // quasi-fermi level, exp(-Ef/Vt)
    complex <double> v;   //                  , exp(Ef/Vt)
    complex <double> r;   // recombination rate, SRH
    complex <double> rho; // charge density, C/cm3
    complex <double> Crho;// charge density, C/cm3
    double mun; // mobility, nm^2/v.s
    double mup;
    double tau;
    complex <double> Ex;  // electric field
    complex <double> Ey;
    complex <double> Ez;
    double Type;// Material Type 1=channel, 2=insulator, 3=electrolyte, 4=sub, 5=analyte;
};

class ACDD : public FVMesh{

public:
    ACDD();
    ~ACDD();

    //Process
    void IdVG2D();
    void IdVD2D();
    void IdVG3D();
    void IdVD3D();
    //void GComparison3D();
    //void GComparison3D_Dotshift();
    //void GComparison3D_Dotshift_2NWR();
    //void RefineIdVG3D(double NewSimTol);
    //void RefineIdVD3D(double NewSimTol);
    //void RefineGComparison3D();
    //void RefineGComparison3D_Dotshift();

    // DD parameter setting
    void ACDDmodelParameterSet();

    //DD Initial Guess
    void DDInitialGuess2D();
    void DDInitialGuessPNJunction2D();
    void DDInitialGuessMOSFET2D();
    void DDInitialGuessISFET2D();

    void DDInitialGuess3D();
    void DDInitialGuessPNJunction3D();
    void DDInitialGuessMOSFET3D();
    void DDInitialGuessISFET3D();
    void DDInitialGuess1NWR3D();
    void DDInitialGuess2NWR3D();

    //Poisson Solver
    double PoissonSolver2D();
    double PoissonGaussSeidel2D();
    complex <double> PoissonGaussSeidelInner2D(int i, int j);

    double PoissonSolver3D();
    double PoissonGaussSeidel3D();
    complex <double> PoissonGaussSeidelInner3D(int i, int j, int k);

    //ElectronContinuity(EC)
    double ECSolver2D();
    double ECTypeA2D();
    double ECTypeB2D();
    complex <double> ECInner2D(int i, int j);

    double ECSolver3D();
    double ECTypeA3D();
    double ECTypeB3D();
    double EC1NWR3D(); // NWR only
    double EC2NWR3D(); // NWR only
    complex <double> ECInner3D(int i, int j, int k);

    //HoleContinuity(HC)
    double HCSolver2D();
    double HCTypeA2D();
    double HCTypeB2D();
    complex <double> HCInner2D(int i, int j);

    double HCSolver3D();
    double HCTypeA3D();
    double HCTypeB3D();
    double HC1NWR3D(); // NWR only
    double HC2NWR3D(); // NWR only
    complex <double> HCInner3D(int i, int j, int k);

    //BC
    void PoissonBC2D();
    void PoissonBC2D_PN();
    void PoissonBC2D_MOSFET();
    void PoissonBC2D_ISFET();
    void ECBC2D();
    void ECBC2D_PN();
    void ECBC2D_MOSFET();
    void ECBC2D_ISFET();
    void HCBC2D();
    void HCBC2D_PN();
    void HCBC2D_MOSFET();
    void HCBC2D_ISFET();

    void PoissonBC3D();
    void PoissonBC3D_PN();
    void PoissonBC3D_MOSFET();
    void PoissonBC3D_ISFET();
    void PoissonBC3D_1NWR();
    void PoissonBC3D_2NWR();
    void ECBC3D();
    void ECBC3D_PN();
    void ECBC3D_MOSFET();
    void ECBC3D_ISFET();
    void ECBC3D_1NWR();
    void ECBC3D_2NWR();
    void HCBC3D();
    void HCBC3D_PN();
    void HCBC3D_MOSFET();
    void HCBC3D_ISFET();
    void HCBC3D_1NWR();
    void HCBC3D_2NWR();

    //Tool Function
    void PrintMaterial2D(string path);
    void PrintMaterial3D(string path);
    void ReadMaterial2D(string path);
    void ReadMaterial3D(string path);

    //Current Calculations
    void Jcal2D ();
    void Jcal2D_PN();
    void Jcal2D_MOSFET();
    void Jcal2D_ISFET();

    void JcalS2D_PN(complex <double> &JSn, complex <double> &JSp);
    void JcalD2D_PN(complex <double> &JDn, complex <double> &JDp);

    void JcalS2D_MOSFET(complex <double> &JSn, complex <double> &JSp);
    void JcalD2D_MOSFET(complex <double> &JDn, complex <double> &JDp);
    void JcalB2D_MOSFET(complex <double> &JBn, complex <double> &JBp);

    void JcalS2D_ISFET(complex <double> &JSn, complex <double> &JSp);
    void JcalD2D_ISFET(complex <double> &JDn, complex <double> &JDp);
    void JcalB2D_ISFET(complex <double> &JBn, complex <double> &JBp);

    void Jcal3D();
    void Jcal3D_PN();
    void Jcal3D_MOSFET();
    void Jcal3D_ISFET();
    void Jcal3D_1NWR();
    void Jcal3D_2NWR();

    void JcalS3D_PN(complex <double> &JSn, complex <double> &JSp);
    void JcalD3D_PN(complex <double> &JDn, complex <double> &JDp);

    void JcalS3D_MOSFET(complex <double> &JSn, complex <double> &JSp);
    void JcalD3D_MOSFET(complex <double> &JDn, complex <double> &JDp);
    void JcalB3D_MOSFET(complex <double> &JBn, complex <double> &JBp);

    void JcalS3D_ISFET(complex <double> &JSn, complex <double> &JSp);
    void JcalD3D_ISFET(complex <double> &JDn, complex <double> &JDp);
    void JcalB3D_ISFET(complex <double> &JBn, complex <double> &JBp);

    void JcalS3D_NWR1(complex <double> &JSn, complex <double> &JSp); //NWR1
    void JcalD3D_NWR1(complex <double> &JDn, complex <double> &JDp);

    void JcalS3D_NWR2(complex <double> &JSn, complex <double> &JSp); //NWR2
    void JcalD3D_NWR2(complex <double> &JDn, complex <double> &JDp);

    //Efield Calculations
    void EfieldCalculation2D();
    void EfieldCalculation3D();

    //Space Charge Calculations
    void RhoCalculation2D();
    void RhoCalculation3D();

    //Parameter Calculation
    double munCal(double T, double dopping, int f);
    double mupCal(double T,double dopping, int f);
    double tauNCal(double dopping);
    double tauPCal(double dopping);
    complex <double> SRHrecomb2D(int i, int j);
    complex <double> SRHrecomb3D(int i, int j, int k);

    //Add Analyte
    void AddDotString3D(double Yshift);
    void AddDot3D(double DotXCenter, double DotYCenter, int AnalyteFlag);

    //Find NWR BC
    void FindNWRBC1();
    void FindNWRBC2();

private:

    //Tool Function
    complex <double> Bern(complex <double>, complex <double>);
    void BernoulliX();
    void DDInitialize();

protected:

    //*********************************
    //Drift Diffusion Model and Poisson
    //*********************************

    Semiconductor_complex *DDmaterial;

    double SimTolPoisson, SimTolEC, SimTolHC;
    int DD_loop, maxIter;
    double bernXl,Na, Nd, Nai, Ndi, SiNc, SiNv, ZnONc, ZnONv;
    double NaPlus, NdPlus, NaPlusi, NdPlusi;
    double Si_ni_m, Si_ni_cm, Si_ni_nm;
    double ZnO_ni_m, ZnO_ni_cm, ZnO_ni_nm;
    double ni_m, ni_cm, ni_nm;
    double mun_semi, mup_semi, mun_sub, mup_sub;
    complex <double> volS, volDe, volDs, volDi, volD, volB, volGe, volGs, volGi, volG, wfG;
    double mun_electrolyte, mup_electrolyte, C0, CC0;
    double AnalyteRadius, AnalyteValence, ShiftDistanceY,ReceptorLength,AnalytePermittivity;
    int SubstrateTOP=0,ElectrolyteBottom=0;
    int NWRleft1=0,NWRright1=0,NWRtop1=0,NWRbottom1=0,NWRcenteryP1=0,NWRcenterzP1=0;
    int NWRleft2=0,NWRright2=0,NWRtop2=0,NWRbottom2=0,NWRcenteryP2=0,NWRcenterzP2=0;
    int DotNumber;

    /*
    int DD_loop,Kicking_loop,SubstrateTOP,ElectrolyteBottom;
    int NWRleft,NWRright,NWRtop,NWRbottom,NWRcenterxP,NWRcenterzP;
    int NWRleft2,NWRright2,NWRtop2,NWRbottom2,NWRcenterxP2,NWRcenterzP2;
    double volS, volDe, volDs, volDi, volD, volB, volGe, volGs, volGi, volG, wfG;
    double sbvn, sbvp;
    double mun_semi, mup_semi, mun_sub, mup_sub;
    double SimTolPoisson, SimTolEC, SimTolHC, SimTolElectron, maxIter;
    double *PsiKick, *ErrorElectronNew, *ErrorElectronPre;
    double *XFinal, *YFinal, *ZFinal;
    double junction_depth, junction_length, tox_thickness;
    double substrate_thickness, box_thickness, channel_thickness, electrolyte_thickness, channel_length;
    double Lambda_D1, Lambda_D2, D_KCl, Omega, CDL, COxide, CTotal, Conductivity,CapLambda;
    double mun_electrolyte, mup_electrolyte, C0, CC0;
    int AnalyteNum;
    */

};

#endif // ACDD_H
