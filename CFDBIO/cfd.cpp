#include <math.h>
#include <complex>
# include <fstream>
# include <iostream>
# include <string>
# include <sstream>

#include "cfd.h"
#include "Parameter.h"


using namespace std;


CFD::CFD(){
}

CFD::~CFD(){
}

void CFD::CFD_NSParameter(){

    /*
     *unit is nm !
     *
     *Dynamic Viscosity (Absolute Viscosity) :  1e-3 Pa s
     *( 1 Pa s = 1 N s/m^2 = 1 kg/(m s) )
     *
     *Water : 1e-3 Pa s  = 1e-3 kg/(m s) = 1e-12 kg/(nm s)
     *
     *Kinematic Viscosity = Dynamic viscosity / rho = m^2/s
     *
     *We are using Dynamic Viscosity here.
     *
     *Water Density 1000 kg/m^3 = 1e-24 kg/nm^3
     *
     *units in this program:
     *Length:    nm
     *Velocity:  nm/s
     *Dendity:   kg/nm^3
     *Viscosity: N*s/nm^2
     *Pressure:  N/nm^2
     */

    SimTolVelocity=1e-5;
    SimTolPressure=1e-30;
    SimTolPoisson=1e-6;

    MaxIter=50000;

    // alpha 0~1  pressure need more test.
    alphadP=1;
    alphaP=1;
    alphaPc=1; //more nodes -> smaller alphaPc

    //Warrning ! 3D SIMPLER will need this Velocity relaxzation!! But 2D won't need!!
    alphaVx=alphaVy=alphaVz=1;

    Pbreak=dPbreak=3;
    Vxbreak=Vybreak=Vzbreak=3;

    Vac=0.5;
    CC0=1;
    D_KCl=2.47e-9; //[m2/s]
    Conductivity=11.824e-3; //[S/m] // 1mM 11.824 mS/m

    //ion concentration
    //10*6.022e-4;   6.022*e-4 1/nm3 = 1mM = 1 mol/m3 = 1e-3 mol/L = 1e-3 M, 3nm debye = 10mM
    //PH=7 => H+=1e-7 M = 1e-4*1e-3 M = 6.022*e-8 1/nm3

    Lambda_D1=sqrt(80*e0*1e9*kb*Tamb/(q0*q0*Avogadro*2*CC0)); //[m]
    //Lambda_D2=sqrt(D_KCl*80*e0*1e9/Conductivity); //[m]
    Lambda_D2=Lambda_D1; //[m]

    //cerr <<"Lambda_D1="<<Lambda_D1*1e9<<" [nm]"<<endl;
    //cerr <<"Lambda_D2="<<Lambda_D2*1e9<<" [nm]"<<endl;

    CDL=Water_permi*e0*1e9/Lambda_D2; //[F/m2]
    COxide=3.9*e0*1e18/4;
    CTotal=CDL*COxide/(CDL+COxide);


    CharFreq=Conductivity*Lambda_D2/(50e-6*Water_permi*e0*1e9);
    //ACFreq=CharFreq*pow(pow(10,0.2),Global_Index_i);
    ACFreq=CharFreq;
    Omega=2*M_PI*ACFreq;
    CapLambda=0.25;

    //cerr <<"ACFreq="<< ACFreq <<"[Hz]"<<endl;
    //cerr <<"BigOCFD="<< OCFD/CharFreq <<endl;
}

void CFD::CFD_NSInitialGuessComplex2D(){

    /*
     *Data points:
     *Total:px
     *Array:[0]-[px-1]
     *Total data points for Pressure:for (int i=0; i<px; i++)
     *Inner data points for Pressure:for (int i=1; i<px-1; i++)
     *Total data points for Velocity:for (int i=1; i<px; i++)
     *Inner data points for Velocity:for (int i=2; i<px-1; i++)
     */

    /*
     *All initial guess region setting should be using <= or >=,
     *and corresponding boundary condition should be using < or >
     */

    u = new Vx [L];
    v = new Vy [L];

    //u boundary guess , i wall
    // Two Electrode
    for (int i=1; i<px; i++) {
        int pointer1,pointer1_in;
        double CoordX, xstep;

        pointer1 = (px)*(0) + (i);
        pointer1_in = (px)*(0) + (i-1);

        CoordX=(mesh[pointer1].coordX+mesh[pointer1_in].coordX)/2;
        xstep=abs(mesh[pointer1].coordX-mesh[pointer1_in].coordX);

        if(CoordX>lx/2+ElectrodeGap/2 && CoordX<lx/2+ElectrodeGap/2+ElectrodeWidth){
            double Vapp=Vac;
            u[pointer1].u=CapLambda*((-1)*80*e0*1e9/(4*eta*1e18)*(norm(CFDmaterial[pointer1].phi-Vapp)-norm(CFDmaterial[pointer1_in].phi-Vapp))/(xstep*1e-9))*1e9;
            //0.25 is capital lambda factoer account for Stern layer
        }

        if(CoordX<lx/2-ElectrodeGap/2 && CoordX>lx/2-ElectrodeGap/2-ElectrodeWidth){
            double Vapp=(-1)*Vac;
            u[pointer1].u=CapLambda*((-1)*80*e0*1e9/(4*eta*1e18)*(norm(CFDmaterial[pointer1].phi-Vapp)-norm(CFDmaterial[pointer1_in].phi-Vapp))/(xstep*1e-9))*1e9;
            //0.25 is capital lambda factoer account for Stern layer
        }
    }
}

void CFD::CFD_SIMPLE2D(){

    int Vxloop(0),Vyloop(0),dPloop(0),NumIter(0);
    double errMax, errVx, errVy, errdP;

    ofstream output2;

    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);
    output2 <<"========================================================"<<endl;

    do{
        errMax=0;

        stringstream nameP1, nameVx1, nameVy1;
        string nameP2, nameVx2, nameVy2;

        NumIter++;

        nameP1<<"SIMPLE_P_"<<NumIter<<".txt";
        nameP2=nameP1.str();

        nameVx1<<"SIMPLE_Vx_"<<NumIter<<".txt";
        nameVx2=nameVx1.str();

        nameVy1<<"SIMPLE_Vy_"<<NumIter<<".txt";
        nameVy2=nameVy1.str();

        //Step 1
        output2 <<"Iter:" << NumIter <<'\t';
        //VxSolver========================
        errVx=CFD_VxSolver2D();
        Vxloop=NSloop;
        if(errVx>errMax)
            errMax=errVx;

        output2 <<"Vx:" << Vxloop <<'\t'<<"errVx:" <<errVx<<'\t';
        if( NumIter==1|| NumIter==2 || NumIter==5 || NumIter==10  || NumIter==20 || NumIter==100 || NumIter==500 || NumIter%1000==0){
            cerr <<"Iter:" << NumIter <<'\t' <<"Vx:" << Vxloop <<'\t'<<"errVx:"<<errVx<<'\t';
        }


        //VySolver========================
        errVy=CFD_VySolver2D();
        Vyloop=NSloop;
        if(errVy>errMax)
            errMax=errVy;

        output2 <<"Vy:" << Vyloop <<'\t'<<"errVy:"<<errVy<<'\t';
        if( NumIter==1|| NumIter==2 || NumIter==5 || NumIter==10  || NumIter==20 || NumIter==100 || NumIter==500 || NumIter%1000==0 ){
            cerr <<"Vy:" << Vyloop <<'\t'<<"errVy:"<<errVy<<'\t';
        }


        //Step 2
        //dPSolver========================
        errdP=CFD_dPSolver2D();
        dPloop=NSloop;
        if(errdP>errMax)
            errMax=errdP;

        output2 <<"dP:" << dPloop <<'\t'<<"errdP:"<<errdP<<'\t';
        output2 <<"errMax:" << errMax <<'\t' <<endl;
        if( NumIter==1|| NumIter==2 || NumIter==5 || NumIter==10  || NumIter==20 || NumIter==100 || NumIter==500 || NumIter%1000==0){
            cerr <<"dP:" << dPloop <<'\t'<<"errdP:"<<errdP<<'\t'<<"errMax:" << errMax <<'\t' <<endl;
        }



        //Step 3
        //Correct Vx Vy dP========================
        CFD_dVxCalculation2D();
        CFD_dVyCalculation2D();

        if(NumIter<=10000){
            if( NumIter==1|| NumIter==2 || NumIter==5  || NumIter==10  || NumIter==20 || NumIter==100 || NumIter==500 || NumIter%1000==0){
                CFD_PrintVx2D(nameVx2.c_str());
                CFD_PrintVy2D(nameVy2.c_str());
                CFD_VelocityCalculation2D();
                CFD_PrintMaterialComplex2D(nameP2.c_str());
            }
        }else{
            if( NumIter%2000==0){
                CFD_PrintVx2D(nameVx2.c_str());
                CFD_PrintVy2D(nameVy2.c_str());
                CFD_VelocityCalculation2D();
                CFD_PrintMaterialComplex2D(nameP2.c_str());
            }
        }

        CFD_VxCorrection2D();
        CFD_VyCorrection2D();
        CFD_PCorrection2D();

    }while( (Vxloop!=1 || Vyloop!=1) && NumIter<MaxIter);


    CFD_VelocityCalculation2D();
    CFD_PrintMaterialComplex2D("SIMPLE_P.dat");
    CFD_PrintVx2D("SIMPLE_Vx.dat");
    CFD_PrintVy2D("SIMPLE_Vy.dat");
}

void CFD::CFD_SIMPLER2D(){

    int Vxloop(0),Vyloop(0),Ploop(0),dPloop(0),NumIter(0);
    double errMax, errVx, errVy, errP, errdP;

    ofstream  output2;

    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);
    output2 <<"========================================================"<<endl;

    do{
        errMax=0;

        stringstream nameP1, nameVx1, nameVy1;
        string nameP2, nameVx2, nameVy2;

        NumIter++;

        nameP1<<"SIMPLER_P_2D_"<<NumIter<<".txt";
        nameP2=nameP1.str();

        nameVx1<<"SIMPLER_Vx_2D_"<<NumIter<<".txt";
        nameVx2=nameVx1.str();

        nameVy1<<"SIMPLER_Vy_2D_"<<NumIter<<".txt";
        nameVy2=nameVy1.str();

        //Step 1
        output2 <<"Iter:" << NumIter <<'\t';
        //Psedo-Vx========================
        CFD_pseudoVx2D();
        //Psedo-Vy========================
        CFD_pseudoVy2D();

        //Step 2
        //PSolver========================
        errP=CFD_PSolver2D();
        Ploop=NSloop;
        if(errP>errMax)
            errMax=errP;

        output2 <<"P:" << Ploop <<'\t'<<errP<<'\t';
        if( NumIter==1|| NumIter==2 || NumIter==5 || NumIter==10  || NumIter==20  || NumIter==50  || NumIter==100  || NumIter==200 || NumIter==500 || NumIter%1000==0 || NumIter%2000==0 ){
            cerr    <<"Iter:" << NumIter <<'\t' <<"P:" << Ploop <<'\t'<<errP<<'\t';
        }

        //Step 3
        //VxSolver========================
        errVx=CFD_VxSolver2D();
        Vxloop=NSloop;
        if(errVx>errMax)
            errMax=errVx;

        output2 <<"Vx:" << Vxloop <<'\t'<<"errVx:" <<errVx<<'\t';
        if( NumIter==1|| NumIter==2 || NumIter==5 || NumIter==10  || NumIter==20  || NumIter==50  || NumIter==100  || NumIter==200 || NumIter==500 || NumIter%1000==0 || NumIter%2000==0 ){
            cerr <<"Vx:" << Vxloop <<'\t'<<"errVx:"<<errVx<<'\t';
        }

        //VySolver========================
        errVy=CFD_VySolver2D();
        Vyloop=NSloop;
        if(errVy>errMax)
            errMax=errVy;

        output2 <<"Vy:" << Vyloop <<'\t'<<"errVy:"<<errVy<<'\t';
        if( NumIter==1|| NumIter==2 || NumIter==5 || NumIter==10  || NumIter==20  || NumIter==50  || NumIter==100  || NumIter==200 || NumIter==500 || NumIter%1000==0 || NumIter%2000==0 ){
            cerr <<"Vy:" << Vyloop <<'\t'<<"errVy:"<<errVy<<'\t';
        }

        //Step 4
        //dPSolver========================
        errdP=CFD_dPSolver2D();
        dPloop=NSloop;
        if(errdP>errMax)
            errMax=errdP;

        output2 <<"dP:" << dPloop <<'\t'<<"errdP:"<<errdP<<'\t';
        output2 <<"errMax:" << errMax <<'\t' <<endl;
        if( NumIter==1|| NumIter==2 || NumIter==5 || NumIter==10  || NumIter==20  || NumIter==50  || NumIter==100  || NumIter==200 || NumIter==500 || NumIter%1000==0 || NumIter%2000==0 ){
            cerr <<"dP:" << dPloop <<'\t'<<"errdP:"<<errdP<<'\t'<<"errMax:" << errMax <<'\t' <<endl;
        }

        //Step 5
        //Correct Vx Vy dP========================
        CFD_dVxCalculation2D();
        CFD_dVyCalculation2D();

        if(NumIter<=10000){
            if( NumIter==1|| NumIter==2 || NumIter==5  || NumIter==10  || NumIter==20 || NumIter==100 || NumIter==500 || NumIter%1000==0){
                CFD_PrintVx2D(nameVx2.c_str());
                CFD_PrintVy2D(nameVy2.c_str());
                CFD_VelocityCalculation2D();
                CFD_PrintMaterialComplex2D(nameP2.c_str());
            }
        }else{
            if( NumIter%2000==0){
                CFD_PrintVx2D(nameVx2.c_str());
                CFD_PrintVy2D(nameVy2.c_str());
                CFD_VelocityCalculation2D();
                CFD_PrintMaterialComplex2D(nameP2.c_str());
            }
        }

        CFD_VxCorrection2D();
        CFD_VyCorrection2D();

    }while( (Vxloop!=1 || Vyloop!=1) && NumIter<MaxIter);


    CFD_VelocityCalculation2D();
    CFD_PrintMaterialComplex2D("SIMPLE_P.dat");
    CFD_PrintVx2D("SIMPLE_Vx.dat");
    CFD_PrintVy2D("SIMPLE_Vy.dat");
}

double CFD::CFD_VxSolver2D(){

    NSloop=0;

    double errVx_new(0),errVx_max(0);

    do{
        NSloop++;

        errVx_new=CFD_VxGaussSeidel2D();

        if(errVx_max < errVx_new) {errVx_max=errVx_new;}

        //if(NSloop==1 || NSloop%500==0)
        //cerr <<"Vx:"<< NSloop <<'\t' << errVx_new<<endl;

        if(NSloop==Vxbreak) break;

    }while(errVx_new>SimTolVelocity);

    return errVx_max;
}

double CFD::CFD_VxGaussSeidel2D(){

    double max_val=0;

#pragma omp parallel for reduction(max : max_val)
    for (int i=2; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double Vx=u[pointer].u;

            u[pointer].u=CFD_VxGaussSeidelInner2D(i,j);

            double error=abs(u[pointer].u-Vx);

            error=error/(abs(Vx)+1);

            if(error > max_val)
                max_val=error;
        }
    }

    CFD_VxBoundary2D();

    return max_val;

}

double CFD::CFD_VxGaussSeidelInner2D(int i, int j){

    int pointer    =   (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_inn=   (px)*(j) + (i-2);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double xstep_p=abs((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2);
    double xstep_n=abs((mesh[pointer].coordX-mesh[pointer_inn].coordX)/2);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer].coordY-mesh[pointer_jn].coordY);

    double deltax=abs(mesh[pointer].coordX-mesh[pointer_in].coordX);
    double deltay=abs((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2);

    double aN=eta*deltax/ystep_p;
    double aS=eta*deltax/ystep_n;
    double aE=eta*deltay/xstep_p;
    double aW=eta*deltay/xstep_n;
    double aP=aN+aS+aE+aW;

    double Sum=aE*u[pointer_ip].u+aW*u[pointer_in].u
              +aN*u[pointer_jp].u+aS*u[pointer_jn].u
              +deltay*(CFDmaterial[pointer_in].p-CFDmaterial[pointer].p);
            //+Source term Bij

    double uk=u[pointer].u;

    return alphaVx*Sum/aP+(1-alphaVx)*uk;
}

void CFD::CFD_dVxCalculation2D(){

#pragma omp parallel for
    for (int i=2; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            int pointer = (px)*(j) + (i);
            int pointer_in =   (px)*(j) + (i-1);
            int pointer_jp =   (px)*(j+1) + (i);
            int pointer_jn =   (px)*(j-1) + (i);

            double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;

            u[pointer].du=deltay/CFD_uAcoef2D(i,j)*(CFDmaterial[pointer_in].dp-CFDmaterial[pointer].dp);

        }
    }

    for (int i=1; i<px; i++) {
        for (int j=0; j<py; j++) {
            int pointer = (px)*(j) + (i);
            u[pointer].cu=u[pointer].u+u[pointer].du;
        }
    }
}

void CFD::CFD_pseudoVx2D(){

    /*
     *Data points:
     *Total:px
     *Array:[0]-[px-1]
     *Total data points for Pressure:for (int i=0; i<px; i++)
     *Inner data points for Pressure:for (int i=1; i<px-1; i++)
     *Total data points for Velocity:for (int i=1; i<px; i++)
     *Inner data points for Velocity:for (int i=2; i<px-1; i++)
     */

    //Copy u to pseudo-u
    for (int i=1; i<px; i++) {
        for (int j=0; j<py; j++) {
            int pointer =(px)*(j) + (i);
            u[pointer].pu=u[pointer].u;
        }
    }

    //Calculate pseudo-u
    for (int i=2; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);
            int pointer_ip =   (px)*(j) + (i+1);
            int pointer_in =   (px)*(j) + (i-1);
            int pointer_inn =  (px)*(j) + (i-2);
            int pointer_jp =   (px)*(j+1) + (i);
            int pointer_jn =   (px)*(j-1) + (i);

            double xstep_p=abs((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2);
            double xstep_n=abs((mesh[pointer].coordX-mesh[pointer_inn].coordX)/2);
            double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
            double ystep_n=abs(mesh[pointer].coordY-mesh[pointer_jn].coordY);
            double deltax=abs(mesh[pointer].coordX-mesh[pointer_in].coordX);
            double deltay=abs((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2);

            double aN=eta*deltax/ystep_p;
            double aS=eta*deltax/ystep_n;
            double aE=eta*deltay/xstep_p;
            double aW=eta*deltay/xstep_n;
            double aP=aN+aS+aE+aW;

            double Sum=aE*u[pointer_ip].pu+aW*u[pointer_in].pu
                      +aN*u[pointer_jp].pu+aS*u[pointer_jn].pu;
                    //+Source term Bij

            u[pointer].pu=Sum/aP;
        }
    }
}

double CFD::CFD_VySolver2D(){

    NSloop=0;

    double errVy_new(0), errVy_max(0);

    do{
        NSloop++;

        errVy_new=CFD_VyGaussSeidel2D();

        if(errVy_max <  errVy_new) {errVy_max= errVy_new;}

        //if(NSloop==1 || NSloop%500==0)
        //cerr <<"Vy:"<< NSloop <<'\t' <<errVy_ratio<<'\t'<<errVy_new<<endl;

        if(NSloop==Vybreak) break;

    }while(errVy_new>SimTolVelocity);

    return errVy_max;

}

double CFD::CFD_VyGaussSeidel2D(){

    double max_val=0;

#pragma omp parallel for reduction(max : max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=2; j<py-1; j++) {

            int pointer_v = (px)*(j) + (i);

            double Vy=v[pointer_v].v;

            v[pointer_v].v=CFD_VyGaussSeidelInner2D(i,j);

            double error=abs(v[pointer_v].v-Vy);

            error=error/(abs(Vy)+1);

            if(error > max_val)
                max_val=error;
        }
    }

    CFD_VyBoundary2D();

    return max_val;

}

double CFD::CFD_VyGaussSeidelInner2D(int i, int j){

    int pointer    =   (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);
    int pointer_jnn=   (px)*(j-2) + (i);

    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer].coordX-mesh[pointer_in].coordX);
    double ystep_p=abs((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2);
    double ystep_n=abs((mesh[pointer].coordY-mesh[pointer_jnn].coordY)/2);
    double deltax=abs((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2);
    double deltay=abs(mesh[pointer].coordY-mesh[pointer_jn].coordY);

    double aN=eta*deltax/ystep_p;
    double aS=eta*deltax/ystep_n;
    double aE=eta*deltay/xstep_p;
    double aW=eta*deltay/xstep_n;
    double aP=aN+aS+aE+aW;

    double Sum=aE*v[pointer_ip].v+aW*v[pointer_in].v
              +aN*v[pointer_jp].v+aS*v[pointer_jn].v
              +deltax*(CFDmaterial[pointer_jn].p-CFDmaterial[pointer].p);
            //+Source term Bij

    double vk=v[pointer].v;

    return alphaVy*Sum/aP+(1-alphaVy)*vk;
}

void CFD::CFD_dVyCalculation2D(){

    /*
     *Data points:
     *Total:px
     *Array:[0]-[px-1]
     *Total data points for Pressure:for (int i=0; i<px; i++)
     *Inner data points for Pressure:for (int i=1; i<px-1; i++)
     *Total data points for Velocity:for (int i=1; i<px; i++)
     *Inner data points for Velocity:for (int i=2; i<px-1; i++)
     */

#pragma omp parallel for
    for (int i=1; i<px-1; i++) {
        for (int j=2; j<py-1; j++) {

            int pointer = (px)*(j) + (i);
            int pointer_ip =   (px)*(j) + (i+1);
            int pointer_in =   (px)*(j) + (i-1);
            int pointer_jn =   (px)*(j-1) + (i);

            double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;

            v[pointer].dv=deltax/CFD_vAcoef2D(i,j)*(CFDmaterial[pointer_jn].dp-CFDmaterial[pointer].dp);

        }
    }
    for (int i=0; i<px; i++) {
        for (int j=1; j<py; j++) {
            int pointer = (px)*(j) + (i);
            v[pointer].cv=v[pointer].v+v[pointer].dv;
        }
    }
}

void CFD::CFD_pseudoVy2D(){

    /*
     *Data points:
     *Total:px
     *Array:[0]-[px-1]
     *Total data points for Pressure:for (int i=0; i<px; i++)
     *Inner data points for Pressure:for (int i=1; i<px-1; i++)
     *Total data points for Velocity:for (int i=1; i<px; i++)
     *Inner data points for Velocity:for (int i=2; i<px-1; i++)
     */

    //Copy v to pseudo-v
    for (int i=0; i<px; i++) {
        for (int j=1; j<py; j++) {
            int pointer =(px)*(j) + (i);
            v[pointer].pv=v[pointer].v;
        }
    }

    //Calculate pseudo-v
    for (int i=1; i<px-1; i++) {
        for (int j=2; j<py-1; j++) {

            int pointer = (px)*(j) + (i);
            int pointer_ip =   (px)*(j) + (i+1);
            int pointer_in =   (px)*(j) + (i-1);
            int pointer_jp =   (px)*(j+1) + (i);
            int pointer_jn =   (px)*(j-1) + (i);
            int pointer_jnn =  (px)*(j-2) + (i);

            double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
            double xstep_n=abs(mesh[pointer].coordX-mesh[pointer_in].coordX);
            double ystep_p=abs((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2);
            double ystep_n=abs((mesh[pointer].coordY-mesh[pointer_jnn].coordY)/2);
            double deltax=abs((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2);
            double deltay=abs(mesh[pointer].coordY-mesh[pointer_jn].coordY);

            double aN=eta*deltax/ystep_p;
            double aS=eta*deltax/ystep_n;
            double aE=eta*deltay/xstep_p;
            double aW=eta*deltay/xstep_n;
            double aP=aN+aS+aE+aW;

            double Sum=aE*v[pointer_ip].v+aW*v[pointer_in].v
                      +aN*v[pointer_jp].v+aS*v[pointer_jn].v;
                    //+Source term Bij

            v[pointer].pv=Sum/aP;
        }
    }
}

double CFD::CFD_dPSolver2D(){

    NSloop=0;

    double errdP_new(0),errdP_max(0);

    do{
        NSloop++;

        errdP_new=CFD_dPGaussSeidel2D();

        if(errdP_max < errdP_new) {errdP_max=errdP_new;}

        //if(NSloop==1 || NSloop%500==0)
        //cerr <<"dP:"<< NSloop <<'\t'<<errdP_new<<endl;

        if(NSloop==dPbreak) break;

    }while(errdP_new>SimTolPressure);

    return errdP_max;

}

double CFD::CFD_dPGaussSeidel2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double dP=CFDmaterial[pointer].dp;

            CFDmaterial[pointer].dp=CFD_dPGaussSeidelInner2D(i,j);

            double error=abs(CFDmaterial[pointer].dp-dP);

            error=error/(abs(dP)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    return max_val;
}

double CFD::CFD_dPGaussSeidelInner2D(int i, int j){

    //int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double deltax=abs((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2);
    double deltay=abs((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2);

    double aW=0,aE=0,aN=0,aS=0,aP=0;

    if(i>1){
        aW=Drho*deltay/CFD_uAcoef2D(i,j)*deltay;   //a(I-1,J)
    }
    if(j>1){
        aS=Drho*deltax/CFD_vAcoef2D(i,j)*deltax;   //a(I,J-1)
    }
    if(i<px-2){
        aE=Drho*deltay/CFD_uAcoef2D(i+1,j)*deltay; //a(I+1,J)
    }
    if(j<py-2){
        aN=Drho*deltax/CFD_vAcoef2D(i,j+1)*deltax; //a(I,J+1)
    }

    aP=aN+aS+aE+aW;

    double Sum= aE*CFDmaterial[pointer_ip].dp+aW*CFDmaterial[pointer_in].dp+
                aN*CFDmaterial[pointer_jp].dp+aS*CFDmaterial[pointer_jn].dp+
                CFD_BPcoef2D(i,j);

    return  Sum/aP;
}

double CFD::CFD_PSolver2D(){

    NSloop=0;

    double errP_new(0),errP_max(0);

    do{
        NSloop++;

        errP_new=CFD_PGaussSeidel2D();

        if(errP_max < errP_new) {errP_max=errP_new;}

        if(NSloop==Pbreak) break;
        //cerr <<"P:"<< NSloop <<'\t' <<errP_new<<endl;

    }while(errP_new>SimTolPressure);

    return errP_max;
}

double CFD::CFD_PGaussSeidel2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double P=CFDmaterial[pointer].p;

            CFDmaterial[pointer].p=CFD_PGaussSeidelInner2D(i,j);

            double error=abs(CFDmaterial[pointer].p-P);

            error=error/(abs(P)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    return max_val;

}

double CFD::CFD_PGaussSeidelInner2D(int i, int j){

    //int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double deltax=abs((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2);
    double deltay=abs((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2);

    double aW=0,aE=0,aN=0,aS=0,aP=0;

    if(i>1){
        aW=Drho*deltay/CFD_uAcoef2D(i,j)*deltay;   //a(I-1,J)
    }
    if(i<px-2){
        aE=Drho*deltay/CFD_uAcoef2D(i+1,j)*deltay; //a(I+1,J)
    }
    if(j>1){
        aS=Drho*deltax/CFD_vAcoef2D(i,j)*deltax;   //a(I,J-1)
    }
    if(j<py-2){
        aN=Drho*deltax/CFD_vAcoef2D(i,j+1)*deltax; //a(I,J+1)
    }

    aP=aN+aS+aE+aW;

    double Sum= aE*CFDmaterial[pointer_ip].p+aW*CFDmaterial[pointer_in].p+
                aN*CFDmaterial[pointer_jp].p+aS*CFDmaterial[pointer_jn].p+
                CFD_Bcoef2D(i,j);

    return Sum/aP;
}



void CFD::CFD_PrintVx2D(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(8);

    output << "X(1)\tY(2)\tVx(3)\tdVx(4)\tpVx(5)\tcVx(6)#"<<endl;
    output << "[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=1;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            int pointer_in =(px)*(j) + (i-1);
            output << (mesh[pointer].coordX+mesh[pointer_in].coordX)/2 << '\t' << mesh[pointer].coordY  << '\t'
                   << u[pointer].u<< '\t'<< u[pointer].du<< '\t'<< u[pointer].pu<< '\t'<< u[pointer].cu<< endl;
        }
    }

    output.close();
}

void CFD::CFD_PrintVy2D(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(8);

    output << "X(1)\tY(2)\tVy(3)\tdVy(4)\tpVy(5)\tcVy(6)#"<<endl;
    output << "[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=1;j<py;j++){

            int pointer =(px)*(j) + (i);
            int pointer_jn =(px)*(j-1) + (i);
            output << mesh[pointer].coordX<< '\t' << (mesh[pointer_jn].coordY+mesh[pointer].coordY)/2<< '\t'
                   << v[pointer].v << '\t'<< v[pointer].dv<< '\t'<< v[pointer].pv<< '\t'<< v[pointer].cv<< endl;
        }
    }

    output.close();
}


void CFD::CFD_VelocityCalculation2D(){
   /*
    *Data points:
    *Total:px
    *Array:[0]-[px-1]
    *Total data points for Pressure:for (int i=0; i<px; i++)
    *Inner data points for Pressure:for (int i=1; i<px-1; i++)
    *Total data points for Velocity:for (int i=1; i<px; i++)
    *Inner data points for Velocity:for (int i=2; i<px-1; i++)
    */

    // u
    for(int i=1;i<px-1;i++){
        for(int j=0;j<py;j++){

            int pointer = (px)*(j) + (i);
            int pointer_ip =   (px)*(j) + (i+1);

            CFDmaterial[pointer].ui=(u[pointer].u+u[pointer_ip].u)/2;
        }
    }

    for(int j=0;j<py;j++){

        int pointer = (px)*(j) + (0);
        int pointer_ip =   (px)*(j) + (1);
        CFDmaterial[pointer].ui=u[pointer_ip].u;

        pointer = (px)*(j) + (px-1);
        CFDmaterial[pointer].ui=u[pointer].u;

    }

    // v
    for(int i=0;i<px;i++){
        for(int j=1;j<py-1;j++){

            int pointer = (px)*(j) + (i);
            int pointer_jp =   (px)*(j+1) + (i);

            CFDmaterial[pointer].vi=(v[pointer].v+v[pointer_jp].v)/2;
        }
    }
    for(int i=0;i<px;i++){

        int pointer = (px)*(0) + (i);
        int pointer_jp =   (px)*(1) + (i);
        CFDmaterial[pointer].vi=v[pointer_jp].v;

        pointer = (px)*(py-1) + (i);
        CFDmaterial[pointer].vi=v[pointer].v;
    }
}

void CFD::CFD_PrintMaterialComplex2D(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tY(2)\tK(3)\tphi_real(4)\tphi_imag(5)\tu(6)\tv(7)\trho(8)\tmun(7)\tmup(10)\tEx(11)\tEy(12)\tp(13)\tdp(14)\tui(15)\tvi(16)\tflag(17)\tCrho(18)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            output << mesh[pointer].coordX << '\t' << mesh[pointer].coordY << '\t'
                   << CFDmaterial[pointer].k << '\t' <<real(CFDmaterial[pointer].phi)<< '\t' <<imag(CFDmaterial[pointer].phi) << '\t'
                   << CFDmaterial[pointer].u << '\t' << CFDmaterial[pointer].v<< '\t'
                   << CFDmaterial[pointer].rho << '\t'<< CFDmaterial[pointer].mun << '\t'<< CFDmaterial[pointer].mup << '\t'
                   << CFDmaterial[pointer].Ex << '\t'<< CFDmaterial[pointer].Ey << '\t'
                   << CFDmaterial[pointer].p << '\t'<< CFDmaterial[pointer].dp << '\t'<< CFDmaterial[pointer].ui << '\t'
                   << CFDmaterial[pointer].vi << '\t'<< CFDmaterial[pointer].flag << '\t'<< CFDmaterial[pointer].Crho <<endl;
        }
    }

    output.close();
}

void CFD::CFD_ReadMaterialComplex2D(const char *path){

    fstream input;

    input.open(path, fstream::in);

    double buffer;
    double R_real;
    double R_imag;


    // skip first 3 line
    input.ignore(256,'#');
    input.ignore(256,'#');
    input.ignore(256,'#');

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            input >> buffer >> buffer >> CFDmaterial[pointer].k ;
            input >> R_real >> R_imag;
            CFDmaterial[pointer].phi=complex<double> (R_real, R_imag);
            input >> CFDmaterial[pointer].u >> CFDmaterial[pointer].v >> CFDmaterial[pointer].rho
                  >> CFDmaterial[pointer].mun >> CFDmaterial[pointer].mup >> CFDmaterial[pointer].Ex >> CFDmaterial[pointer].Ey
                  >> CFDmaterial[pointer].p >> CFDmaterial[pointer].dp >> CFDmaterial[pointer].ui
                  >> CFDmaterial[pointer].vi >> CFDmaterial[pointer].flag >> CFDmaterial[pointer].Crho;

        }
    }

    input.close();
}

void CFD::CFD_VxCorrection2D(){

    for (int i=2; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            int pointer = (px)*(j) + (i);
            u[pointer].u=u[pointer].cu;
        }
    }
}

void CFD::CFD_VyCorrection2D(){
    for (int i=1; i<px-1; i++) {
        for (int j=2; j<py-1; j++) {
            int pointer = (px)*(j) + (i);
            v[pointer].v=v[pointer].cv;
        }
    }
}

void CFD::CFD_PCorrection2D(){
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            int pointer = (px)*(j) + (i);
            CFDmaterial[pointer].p=CFDmaterial[pointer].p+alphaPc*CFDmaterial[pointer].dp;
        }
    }
}

void CFD::CFD_VxBoundary2D(){ // all boundary conditions are not finished yet

    if(ACDCPhaseFlag==1){
        CFD_VxBoundaryDC2D();
    }
    if(ACDCPhaseFlag==2){
        CFD_VxBoundaryAC2D();
    }
}

void CFD::CFD_VyBoundary2D(){ // all boundary conditions are not finished yet

    if(ACDCPhaseFlag==1){
        CFD_VyBoundaryDC2D();
    }
    if(ACDCPhaseFlag==2){
        CFD_VyBoundaryAC2D();
    }
}

void CFD::CFD_VxBoundaryDC2D(){ // all boundary conditions are not finished yet

}

void CFD::CFD_VyBoundaryDC2D(){ // all boundary conditions are not finished yet

}

void CFD::CFD_VxBoundaryAC2D(){ // all boundary conditions are not finished yet

}

void CFD::CFD_VyBoundaryAC2D(){ // all boundary conditions are not finished yet

}

double CFD::CFD_BPcoef2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;

    return  Drho*(u[pointer].u*deltay-u[pointer_ip].u*deltay+
                  v[pointer].v*deltax-v[pointer_jp].v*deltax);
}

double CFD::CFD_Bcoef2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;


    return  Drho*(u[pointer].pu*deltay-u[pointer_ip].pu*deltay+
                  v[pointer].pv*deltax-v[pointer_jp].pv*deltax);
}

double CFD::CFD_uAcoef2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_inn =  (px)*(j) + (i-2);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double xstep_p=abs((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2);
    double xstep_n=abs((mesh[pointer].coordX-mesh[pointer_inn].coordX)/2);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer].coordY-mesh[pointer_jn].coordY);
    double deltax=abs(mesh[pointer].coordX-mesh[pointer_in].coordX);
    double deltay=abs((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2);

    double aN=eta*deltax/ystep_p;
    double aS=eta*deltax/ystep_n;
    double aE=eta*deltay/xstep_p;
    double aW=eta*deltay/xstep_n;
    double aP=aN+aS+aE+aW;

    return aP;
}

double CFD::CFD_vAcoef2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);
    int pointer_jnn =  (px)*(j-2) + (i);

    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer].coordX-mesh[pointer_in].coordX);
    double ystep_p=abs((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2);
    double ystep_n=abs((mesh[pointer].coordY-mesh[pointer_jnn].coordY)/2);
    double deltax=abs((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2);
    double deltay=abs(mesh[pointer].coordY-mesh[pointer_jn].coordY);

    double aN=eta*deltax/ystep_p;
    double aS=eta*deltax/ystep_n;
    double aE=eta*deltay/xstep_p;
    double aW=eta*deltay/xstep_n;
    double aP=aN+aS+aE+aW;

    return aP;
}

void CFD::CFD_Initializeation2D(){
    for(int i=0;i<L;i++){
        CFDmaterial[i].k=80;
    }
}

void CFD::CFD_ACPoissonInitialGuess2D(){

    CFDmaterial = new CFDElectrolyte [L];

    CFD_Initializeation2D();

    //2D  Two Electrode
    for(int i=0;i<px;i++){
        int pointer1 =(px)*(0) + (i);
        int pointer2 =(px)*(1) + (i);
        //Calcuation is in unit [m]
        if(mesh[pointer1].coordX>=lx/2+ElectrodeGap/2 && mesh[pointer1].coordX<=lx/2+ElectrodeGap/2+ElectrodeWidth){
            double ystep=abs(mesh[pointer1].coordY-mesh[pointer2].coordY)*1e-9;
            double Vapp=Vac;
            complex<double> I(0,1);
            complex<double> A=Conductivity*CFDmaterial[pointer2].phi+Omega*CDL*ystep*Vapp*I;
            complex<double> B(Conductivity,Omega*CDL*ystep);
            CFDmaterial[pointer1].phi=A/B;
            //CFDmaterial[pointer1].phi=Vapp+I*Omega*CDL*Conductivity*(CFDmaterial[pointer2].phi-CFDmaterial[pointer1].phi)/ystep;
        }

        if(mesh[pointer1].coordX<=lx/2-ElectrodeGap/2 && mesh[pointer1].coordX>=lx/2-ElectrodeGap/2-ElectrodeWidth){
            double ystep=abs(mesh[pointer1].coordY-mesh[pointer2].coordY)*1e-9;
            double Vapp=(-1)*Vac;
            complex<double> I(0,1);
            complex<double> A=Conductivity*CFDmaterial[pointer2].phi+Omega*CDL*ystep*Vapp*I;
            complex<double> B(Conductivity,Omega*CDL*ystep);
            CFDmaterial[pointer1].phi=A/B;
            //CFDmaterial[pointer1].phi=Vapp+I*Omega*CDL*Conductivity*(CFDmaterial[pointer2].phi-CFDmaterial[pointer1].phi)/ystep;
        }

        int pointer3 =(px)*(py-1) + (i);
        if(mesh[pointer3].coordX<=lx/2+ReferenceGateWidth/2 && mesh[pointer3].coordX>=lx/2-ReferenceGateWidth/2){
            CFDmaterial[pointer3].phi=0;
        }
    }
}

double CFD::CFD_PoissonSolverComplex2D(){

    int loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        loop++;

        errPhi=CFD_PoissonGaussSeidelComplex2D();

        if(errPhi_max < errPhi) {errPhi_max=errPhi;}

        if(loop<20 || loop%100==0)
        cerr <<"PS:"<< loop <<"\t" <<errPhi<<"\t"<<errPhi_max<<endl;


        if(loop%10000==0)
            CFD_PrintMaterialComplex2D("Complex_temp.txt");

        if(loop==20000)
            break;

    }while(errPhi>SimTolPoisson);

    return errPhi_max;

}

double CFD::CFD_PoissonGaussSeidelComplex2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            complex<double> phik=CFDmaterial[pointer].phi;

            CFDmaterial[pointer].phi=CFD_PoissonGaussSeidelInnerComplex2D(i,j);

            double error=abs(CFDmaterial[pointer].phi-phik);

            error=error/(abs(phik)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    ACPoissonBC2D();

    return max_val;

}

complex<double> CFD::CFD_PoissonGaussSeidelInnerComplex2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double permitivity_ip=CFDmaterial[pointer].k*CFDmaterial[pointer_ip].k / (0.5*CFDmaterial[pointer_ip].k+0.5*CFDmaterial[pointer].k);
    double permitivity_in=CFDmaterial[pointer].k*CFDmaterial[pointer_in].k / (0.5*CFDmaterial[pointer_in].k+0.5*CFDmaterial[pointer].k);
    double permitivity_jp=CFDmaterial[pointer].k*CFDmaterial[pointer_jp].k / (0.5*CFDmaterial[pointer_jp].k+0.5*CFDmaterial[pointer].k);
    double permitivity_jn=CFDmaterial[pointer].k*CFDmaterial[pointer_jn].k / (0.5*CFDmaterial[pointer_jn].k+0.5*CFDmaterial[pointer].k);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);

    double f,df;
    //double volume=deltax*deltay;

    complex<double> phik=CFDmaterial[pointer].phi;

    f=0;
    df=0;

    return (((permitivity_ip*CFDmaterial[pointer_ip].phi/xstep_p+permitivity_in*CFDmaterial[pointer_in].phi/xstep_n)*deltay
            +(permitivity_jp*CFDmaterial[pointer_jp].phi/ystep_p+permitivity_jn*CFDmaterial[pointer_jn].phi/ystep_n)*deltax
            + f - df*phik )
            /
            ((permitivity_ip/xstep_p+permitivity_in/xstep_n)*deltay
            +(permitivity_jp/ystep_p+permitivity_jn/ystep_n)*deltax - df));
}

void CFD::ACPoissonBC2D(){

    //AC Two Electrode
    for(int i=0;i<px;i++){
        int pointer1 =(px)*(0) + (i);
        int pointer2 =(px)*(1) + (i);

        CFDmaterial[pointer1].phi=CFDmaterial[pointer2].phi;

        //Calcuation is in unit [m]
        if(mesh[pointer1].coordX>=lx/2+ElectrodeGap/2 && mesh[pointer1].coordX<=lx/2+ElectrodeGap/2+ElectrodeWidth){
            double ystep=abs(mesh[pointer1].coordY-mesh[pointer2].coordY)*1e-9;
            double Vapp=Vac;
            complex<double> I(0,1);
            complex<double> A=Conductivity*CFDmaterial[pointer2].phi+Omega*CDL*ystep*Vapp*I;
            complex<double> B(Conductivity,Omega*CDL*ystep);
            CFDmaterial[pointer1].phi=A/B;
        }

        if(mesh[pointer1].coordX<=lx/2-ElectrodeGap/2 && mesh[pointer1].coordX>=lx/2-ElectrodeGap/2-ElectrodeWidth){
            double ystep=abs(mesh[pointer1].coordY-mesh[pointer2].coordY)*1e-9;
            double Vapp=(-1)*Vac;
            complex<double> I(0,1);
            complex<double> A=Conductivity*CFDmaterial[pointer2].phi+Omega*CDL*ystep*Vapp*I;
            complex<double> B(Conductivity,Omega*CDL*ystep);
            CFDmaterial[pointer1].phi=A/B;
        }


        if(mesh[pointer1].coordX<=lx/2-ReferenceGateWidth/2 || mesh[pointer1].coordX>=lx/2+ReferenceGateWidth/2){
            pointer1 =(px)*(py-1) + (i);
            pointer2 =(px)*(py-2) + (i);
            CFDmaterial[pointer1].phi=CFDmaterial[pointer2].phi;
        }
    }

    for(int j=0;j<py;j++){
        int pointer1 =(px)*(j) + (0);
        int pointer2 =(px)*(j) + (1);
        CFDmaterial[pointer1].phi=CFDmaterial[pointer2].phi;

        pointer1 =(px)*(j) + (px-1);
        pointer2 =(px)*(j) + (px-2);
        CFDmaterial[pointer1].phi=CFDmaterial[pointer2].phi;
    }
}









