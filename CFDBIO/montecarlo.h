#ifndef MONTECARLO_H
#define MONTECARLO_H
#include "cfd.h"
#include <complex>

struct receptor {
    double coordX;
    double coordY;
    double A;      //concentration
    double B;      //concentration
    double AB;      //concentration
    int flag;
};

struct target {
    double coordX;
    double coordY;
    double coordZ;
    int flag;
};

class MonteCarlo : public CFD{

public:

    MonteCarlo();
    ~MonteCarlo();

    void MC_ParticleTracingParameter3D();
    void MC_ParticleTracingNew();
    void MC_ParticleTracingInitialize3D();
    void MC_ParticleTracingSimulation3D();
    void MC_DirectGenerateOnSurface();
    void MC_Distribution (void);

private:
    //Permutation
    double MC_factorial(double n);
    double MC_Ccomb(double m, double n);
    double MC_Pbinominal(double m, double n, double p);
    double MC_Hcomb(double m, double n);

    //Random Number
    double MC_NormalDistribution();

    //Find IJ for First
    int MC_BCFindI(int dn, int tn);
    int MC_BCFindJ(int dn, int tn);
    int MC_ParticleTracing3DFindI(int dn, int tn);
    int MC_ParticleTracing3DFindJ(int dn, int tn);
    int MC_ParticleTracing3DFindK(int dn, int tn);
    int MC_ParticleTracing3DFindJ_In2DField(int dn, int tn);

    //Capture release function
    void MC_FirstOrderBoundary(int i, int j);
    void MC_FirstOrderRelease(int i, int j);
    void MC_FirstOrderBoundary_FiniteR(int i, int j);
    void MC_FirstOrderRelease_FiniteR(int i, int j);

    //Capture release area
    void MC_StickBoundary_onSensor(int i, int j);
    void MC_RSABoundary_onSensor(int i, int j);
    void MC_FirstOrderBoundary_RSA_onSensor_BlockReceptor(int i, int j);
    void MC_FirstOrderBoundary_RSA_onSensor_Cylindrical(int i, int j);

    void MC_StickBoundary_BottomSurface(int i, int j);
    void MC_RSABoundary_BottomSurface(int i, int j);

    //I/O
    void MC_PrintTrojactory(int tail);


protected:

    //*********************************
    //Brownian & Particle Tracing
    //*********************************

    receptor *ReceptorArray;
    target **Analyte;
    int *ParticleOnWire;

    int TN, RL, DotN, Rpx, Rpy, Mode,WireN,ModeT;
    int MaximumCaptureNumberRSA2;
    int OutputTime, OutputMoleculeNumber, OutputTimeStepsPerFrame, ReportTimeSection;
    double t, tau, AnalyteRadius_m, AnalyteRadius_nm, AnalyteMass, AnalyteValence;
    double D, gamma, tauB, SigmaDist, SigmaForce, g, sigma, mean, Wb;
    double SensorWidth, SensorLength, ReceptorMesh;
    double WireWidth, WireSpacing;
    double C_A, C_B, C_AB;
    double k_forward,k_backward,k_forward1,k_backward1,k_forward2,k_backward2;
    double UnitCon,LC_A,LC_A_cylindrical,R1N_cylindrical;
    double ZetaV;
    /*
    int TN, DotN, RN, AnalyteValence, *TotalCaptureNumber, Rpx, Rpy;
    double UbX, LbX, UbY, LbY, UbZ, LbZ, MaxSpeedx, MaxSpeedy, MaxSpeedz,MaxStepy,MinStepy,MaxStepx,MinStepx,MaxStepz,MinStepz;
    double t ,tau , sigma, mean, D, D_nm, simD, eta_m, AnalyteRadius_m, kd, E, mass, f, g, gamma, tauB, Wb, AnalyteRadius_nm, ReceptorDistance;
    double ReceptorRadius;
    double *Vxi, *Vyi, *Vzi;
    double **DispX, **DispY, **DispZ;
    double **ForceX, **ForceY, **ForceZ;
    double **AccelerationX, **AccelerationY, **AccelerationZ;
    double **VelocityX, **VelocityY, **VelocityZ;
    double **DsquaredX, **DsquaredY, **DsquaredZ;
    double *DsquaredAvg;
    double *AvgStatistic;
    */


};

#endif // MONTECARLO_H
