#ifndef DDMODEL_H
#define DDMODEL_H
#include "fvmesh.h"

struct Semiconductor {
    double k;   // permitivitty
    double dop; // doping, -Nai Ndi
    double phi; // potential, Ei
    double u;   // quasi-fermi level, exp(-Ef/Vt)
    double v;   //                  , exp(Ef/Vt)
    double r;   // recombination rate, SRH
    double rho; // charge density, C/cm3
    double Crho;// charge density, C/cm3
    double mun; // mobility, nm^2/v.s
    double mup;
    double tau;
    double Ex;  // electric field
    double Ey;
    double Ez;
    double Type;// Material Type 1=channel, 2=insulator, 3=electrolyte, 4=sub, 5=analyte;
};

class DDmodel : public FVMesh{

public:
    DDmodel();
    ~DDmodel();

    //Process
    void DD_IdVG2D();
    void DD_IdVD2D();
    void DD_IdVG3D();
    void DD_IdVD3D();
    //void GComparison3D();
    //void GComparison3D_Dotshift();
    //void GComparison3D_Dotshift_2NWR();
    //void RefineIdVG3D(double NewSimTol);
    //void RefineIdVD3D(double NewSimTol);
    //void RefineGComparison3D();
    //void RefineGComparison3D_Dotshift();

    // DD parameter setting
    void DD_ParameterSet();

    //DD Initial Guess
    void DD_NewAndInitialize();
    void DD_InitialGuess2D();
    void DD_InitialGuessPNJunction2D();
    void DD_InitialGuessMOSFET2D();
    void DD_InitialGuessISFET2D();

    void DD_InitialGuess3D();
    void DD_InitialGuessPNJunction3D();
    void DD_InitialGuessMOSFET3D();
    void DD_InitialGuessISFET3D();
    void DD_InitialGuess1NWR3D();
    void DD_InitialGuess2NWR3D();

    //Poisson Solver
    double DD_PoissonSolver2D();
    double DD_PoissonSolver3D();

    //ElectronContinuity(EC)
    double DD_ECSolver2D();
    double DD_ECSolver3D();

    //HoleContinuity(HC)
    double DD_HCSolver2D();
    double DD_HCSolver3D();

    //Tool Function
    void DD_PrintMaterial2D(string path);
    void DD_PrintMaterial3D(string path);
    void DD_ReadMaterial2D(string path);
    void DD_ReadMaterial3D(string path);

    //Add Analyte
    void DD_AddDotString3D(double Yshift);
    void DD_AddDot3D(double DotXCenter, double DotYCenter, int AnalyteFlag);

private:
    //Poisson Solver
    double DD_PoissonGaussSeidel2D();
    double DD_PoissonGaussSeidelInner2D(int i, int j);
    double DD_PoissonGaussSeidel3D();
    double DD_PoissonGaussSeidelInner3D(int i, int j, int k);

    //ElectronContinuity(EC)
    double DD_ECTypeA2D();
    double DD_ECTypeB2D();
    double DD_ECInner2D(int i, int j);
    double DD_ECTypeA3D();
    double DD_ECTypeB3D();
    double DD_EC1NWR3D(); // NWR only
    double DD_EC2NWR3D(); // NWR only
    double DD_ECInner3D(int i, int j, int k);

    //HoleContinuity(HC)
    double DD_HCTypeA2D();
    double DD_HCTypeB2D();
    double DD_HCInner2D(int i, int j);
    double DD_HCTypeA3D();
    double DD_HCTypeB3D();
    double DD_HC1NWR3D(); // NWR only
    double DD_HC2NWR3D(); // NWR only
    double DD_HCInner3D(int i, int j, int k);

    //BC
    void DD_PoissonBC2D();
    void DD_PoissonBC2D_PN();
    void DD_PoissonBC2D_MOSFET();
    void DD_PoissonBC2D_ISFET();
    void DD_ECBC2D();
    void DD_ECBC2D_PN();
    void DD_ECBC2D_MOSFET();
    void DD_ECBC2D_ISFET();
    void DD_HCBC2D();
    void DD_HCBC2D_PN();
    void DD_HCBC2D_MOSFET();
    void DD_HCBC2D_ISFET();
    void DD_PoissonBC3D();
    void DD_PoissonBC3D_PN();
    void DD_PoissonBC3D_MOSFET();
    void DD_PoissonBC3D_ISFET();
    void DD_PoissonBC3D_1NWR();
    void DD_PoissonBC3D_2NWR();
    void DD_ECBC3D();
    void DD_ECBC3D_PN();
    void DD_ECBC3D_MOSFET();
    void DD_ECBC3D_ISFET();
    void DD_ECBC3D_1NWR();
    void DD_ECBC3D_2NWR();
    void DD_HCBC3D();
    void DD_HCBC3D_PN();
    void DD_HCBC3D_MOSFET();
    void DD_HCBC3D_ISFET();
    void DD_HCBC3D_1NWR();
    void DD_HCBC3D_2NWR();

    //Current Calculations
    void DD_Jcal2D ();
    void DD_Jcal2D_PN();
    void DD_Jcal2D_MOSFET();
    void DD_Jcal2D_ISFET();

    void DD_JcalS2D_PN(double &JSn, double &JSp);
    void DD_JcalD2D_PN(double &JDn, double &JDp);

    void DD_JcalS2D_MOSFET(double &JSn, double &JSp);
    void DD_JcalD2D_MOSFET(double &JDn, double &JDp);
    void DD_JcalB2D_MOSFET(double &JBn, double &JBp);

    void DD_JcalS2D_ISFET(double &JSn, double &JSp);
    void DD_JcalD2D_ISFET(double &JDn, double &JDp);
    void DD_JcalB2D_ISFET(double &JBn, double &JBp);

    void DD_Jcal3D();
    void DD_Jcal3D_PN();
    void DD_Jcal3D_MOSFET();
    void DD_Jcal3D_ISFET();
    void DD_Jcal3D_1NWR();
    void DD_Jcal3D_2NWR();

    void DD_JcalS3D_PN(double &JSn, double &JSp);
    void DD_JcalD3D_PN(double &JDn, double &JDp);

    void DD_JcalS3D_MOSFET(double &JSn, double &JSp);
    void DD_JcalD3D_MOSFET(double &JDn, double &JDp);
    void DD_JcalB3D_MOSFET(double &JBn, double &JBp);

    void DD_JcalS3D_ISFET(double &JSn, double &JSp);
    void DD_JcalD3D_ISFET(double &JDn, double &JDp);
    void DD_JcalB3D_ISFET(double &JBn, double &JBp);

    void DD_JcalS3D_NWR1(double &JSn, double &JSp); //NWR1
    void DD_JcalD3D_NWR1(double &JDn, double &JDp);

    void DD_JcalS3D_NWR2(double &JSn, double &JSp); //NWR2
    void DD_JcalD3D_NWR2(double &JDn, double &JDp);

    //Efield Calculations
    void DD_EfieldCalculation2D();
    void DD_EfieldCalculation3D();

    //Space Charge Calculations
    void DD_RhoCalculation2D();
    void DD_RhoCalculation3D();

    //Parameter Calculation
    double DD_munCal(double T, double dopping, int f);
    double DD_mupCal(double T,double dopping, int f);
    double DD_tauNCal(double dopping);
    double DD_tauPCal(double dopping);
    double DD_SRHrecomb2D(int i, int j);
    double DD_SRHrecomb3D(int i, int j, int k);

    //Find NWR BC
    void DD_FindNWRBC1();
    void DD_FindNWRBC2();

    //Tool Function
    double DD_Bern(double, double);
    void DD_BernoulliX();
    void DD_Initialize();

protected:

    //*********************************
    //Drift Diffusion Model and Poisson
    //*********************************

    Semiconductor *DDmaterial;

    double SimTolPoisson, SimTolEC, SimTolHC;
    int DD_loop, maxIter;
    double bernXl,Na, Nd, Nai, Ndi, SiNc, SiNv, ZnONc, ZnONv;
    double NaPlus, NdPlus, NaPlusi, NdPlusi;
    double Si_ni_m, Si_ni_cm, Si_ni_nm;
    double ZnO_ni_m, ZnO_ni_cm, ZnO_ni_nm;
    double ni_m, ni_cm, ni_nm;
    double mun_semi, mup_semi, mun_sub, mup_sub;
    double volS, volDe, volDs, volDi, volD, volB, volGe, volGs, volGi, volG, wfG;
    double mun_electrolyte, mup_electrolyte, C0, CC0;
    double AnalyteRadius, AnalyteValence, ShiftDistanceY,ReceptorLength,AnalytePermittivity;
    int SubstrateTOP=0,ElectrolyteBottom=0;
    int NWRleft1=0,NWRright1=0,NWRtop1=0,NWRbottom1=0,NWRcenteryP1=0,NWRcenterzP1=0;
    int NWRleft2=0,NWRright2=0,NWRtop2=0,NWRbottom2=0,NWRcenteryP2=0,NWRcenterzP2=0;
    int DotNumber;

};


#endif // DDMODEL_H
