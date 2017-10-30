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

    int MC_ParticleTracing3DFindI(int dn, int tn);
    int MC_ParticleTracing3DFindJ(int dn, int tn);
    int MC_ParticleTracing3DFindK(int dn, int tn);

    double MC_NormalDistribution();

    void MC_StickBoundary(int i, int j);
    void MC_RSABoundary(int i, int j);
    void MC_FirstOrderBoundary(int i, int j);
    void MC_FirstOrderBoundary_FiniteR(int i, int j);
    void MC_FirstOrderBoundary_RSA1(int i, int j);
    void MC_FirstOrderBoundary_RSA2(int i, int j);
    int MC_FirstOrderFindI(int dn, int tn);
    int MC_FirstOrderFindJ(int dn, int tn);
    void MC_FirstOrderRelease(int i, int j);
    void MC_FirstOrderRelease_FiniteR(int i, int j);

    double MC_factorial(double n);
    double MC_Ccomb(double m, double n);
    double MC_Pbinominal(double m, double n, double p);
    double MC_Hcomb(double m, double n);
    void MC_Distribution (void);


protected:

    //*********************************
    //Brownian & Particle Tracing
    //*********************************

    receptor *ReceptorArray;
    target **Analyte;
    int *ParticleOnWire;

    int TN, RL, DotN, Rpx, Rpy, ModeT, OutputTimeStep, OutputMoleculeNumber,WireN, OutputFrame;
    int MaximumCaptureNumberRSA2;
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
