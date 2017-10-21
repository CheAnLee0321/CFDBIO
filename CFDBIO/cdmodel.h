#ifndef CDMODEL_H
#define CDMODEL_H
#include "fvmesh.h"

struct electrolyte {
    double k;   // permitivitty
    double C_new;
    double C_old;
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

    void CDInitialize();
    //Convection Diffusion Eq.
    void CDParameter();
    void CDAddReceptors2D();
    void CDInitialGuess2D();
    void CDSolver2D();
    void CDinner2D(int i, int j);
    void CDBoundary2D(int i, int j);

    //Tool Function
    void PrintMaterial2D(string path);
    void PrintReceptor2D(string path);

protected:

    Receptor *ReceptorArray;
    Electrolyte *material;

    int TN, DimensionTag, L1,frame,frame_step;
    double ElectrodeGap, ElectrodeWidth;
    double t, tau, AnalyteRadius_m, AnalyteRadius_nm, eta_m, D,C0;
    double k_forward, k_backward,k_forward1, k_backward1,k_forward2, k_backward2;
    double D_nm;


};

#endif // CDMODEL_H
