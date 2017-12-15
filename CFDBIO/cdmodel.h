#ifndef CDMODEL_H
#define CDMODEL_H
#include "fvmesh.h"
#include "cfd.h"
#include <complex>

struct CDElectrolyte {
    double C_new;
    double C_old;
    double flag;// 1=channel, 2=insulator, 3=electrolyte, 4=sub, 5=analyte;
};

struct Receptor {
    double coordX;
    double coordY;
    double B;      //concentration
    double AB;      //concentration
    int flag;
};

class CDmodel : public CFD{

public:
    CDmodel();
    ~CDmodel();

    void CD_Initialize();
    //Convection Diffusion Eq.
    void CD_Parameter();
    void CD_AddReceptors2D();
    void CD_InitialGuess2D();
    void CD_Solver2D();
    void CD_inner2D(int i, int j);
    void CD_Boundary2D(int i, int j);

    //Tool Function
    void CD_PrintMaterial2D(string path);
    void CD_PrintReceptor2D(string path);

protected:

    Receptor *ReceptorArray;
    CDElectrolyte *CDmaterial;

    int TN, DimensionTag, L1,frame,frame_step;
    double ElectrodeGap, ElectrodeWidth;
    double t, tau, AnalyteRadius_m, AnalyteRadius_nm, eta_m, D,C0;
    double k_forward, k_backward,k_forward1, k_backward1,k_forward2, k_backward2;
    double D_nm;


};

#endif // CDMODEL_H
