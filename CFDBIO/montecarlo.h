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

    void ParticleTracingParameter3D();
    void ParticleTracingNew();
    void ParticleTracingInitialize3D();
    void ParticleTracingSimulation3D();

    void DirectGenerateOnSurface();

    int ParticleTracing3DFindI(int dn, int tn);
    int ParticleTracing3DFindJ(int dn, int tn);
    int ParticleTracing3DFindK(int dn, int tn);

    double NormalDistribution();

    void StickBoundary(int i, int j);
    void RSABoundary(int i, int j);
    void FirstOrderBoundary(int i, int j);
    void FirstOrderBoundary_FiniteR(int i, int j);
    void FirstOrderBoundary_RSA1(int i, int j);
    void FirstOrderBoundary_RSA2(int i, int j);
    int FirstOrderFindI(int dn, int tn);
    int FirstOrderFindJ(int dn, int tn);
    void FirstOrderRelease(int i, int j);
    void FirstOrderRelease_FiniteR(int i, int j);

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
