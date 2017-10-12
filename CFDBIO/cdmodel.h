#ifndef CDMODEL_H
#define CDMODEL_H
#include "fvmesh.h"

struct electrolyte {
    double coordX; //coordinate
    double coordY;
    double coordZ;
    double k;   // permitivitty
    double dop; // doping, -Nai Ndi
    double phi; // potential, Ei
    double u;   // quasi-fermi level, exp(-Ef/Vt)
    double v;   //                  , exp(Ef/Vt)
    double r;   // recombination rate, SRH
    double rho; // charge density, C/cm3
    double Crho; // charge density, C/cm3
    double mun; // mobility, nm^2/v.s
    double mup;
    double tau;
    double Ex;  // electric field
    double Ey;
    double Ez;
    double ui;  //velocity
    double vi;
    double wi;
    double flag;// 1=channel, 2=insulator, 3=electrolyte, 4=sub, 5=analyte;
};

struct receptor {
    double coordX;
    double coordY;
    double B;      //concentration
    double AB;      //concentration
    int flag;
};

typedef receptor Receptor;
typedef electrolyte Electrolyte;

class CDmodel : public FVMesh{

public:
    CDmodel();
    ~CDmodel();

    //Convection Diffusion Eq.
    void CDParameter();
    void CDParameterSet2D();
    void CDStructureParameter2D();
    void CDAddReceptors2D();
    void CDInitialGuess2D();
    void CDSolver2D();
    void CDinner2D(int i, int j);
    void CDBoundary2D(int i, int j);

    void CDParameterSet3D();
    void CDStructureParameter3D();
    void CDAddReceptors3D();
    void CDInitialGuess3D();
    void CDSolver3D();
    void CDinner3D(int i, int j);

protected:

    Receptor *ReceptorArray;
    Electrolyte *material;

    int TN, DimensionTag, L1;
    double ElectrodeGap, ElectrodeWidth;;
    double t, tau, AnalyteRadius_m, AnalyteRadius_nm, eta_m, D;
    double k_forward, k_backward,k_forward1, k_backward1,k_forward2, k_backward2;
    double D_nm;


};

#endif // CDMODEL_H
