#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include "cdmodel.h"
#include "Parameter.h"


CDmodel::CDmodel()
{

}

CDmodel::~CDmodel()
{

}

void CDmodel::CD_Parameter(){

    t=50;       //time, s.
    tau=1e-2;  //time step0
    TN=t/tau+1;        //nubmer of steps
    frame=100;
    frame_step=t/tau/100;

    C0=1;

    AnalyteRadius_m=5e-9;         //radius [m]
    AnalyteRadius_nm=AnalyteRadius_m*1e9;

    eta_m=eta*1e18; //1e-3 [N*s/m^2]
    D=kb*Tamb/(6*M_PI*eta_m*AnalyteRadius_m);
    D_nm=D*1e18;

    k_forward=0;
    k_backward=0;

    k_forward1=1e3;
    k_backward1=1e-3;

    k_forward2=0;
    k_backward2=1e4;


    fstream output;
    output.open("CDParameter2D.txt", fstream::out | fstream::trunc);
    output << "lx="<<lx<< " ly="<<ly<< " lz="<<lz<< endl;
    output << "t="<<t<< " tau="<<tau<<endl;
    output << "AnalyteRadius_m="<<AnalyteRadius_m<<endl;
    output << "Diffusion Constant="<<D<<" (m2/s) "<<D_nm<<" (nm2/s) "<<endl;
    output << "k_forward1="<<k_forward1 << " k_backward1="<<k_backward1<<endl;
    output << "k_forward2="<<k_forward2 << " k_backward2="<<k_backward2<<endl;

    output << "Mx="<<Mx<< " My="<<My<< endl;

    for(int i=0;i<Mx+1;i++){
        output <<"xpin["<<i<<"]="<< xpin[i]<<" ";
    }
    output <<endl;
    for(int i=0;i<My+1;i++){
        output <<"ypin["<<i<<"]="<< ypin[i]<<" ";
    }
    output <<endl;


    for(int i=0;i<Mx;i++){
        output <<"meshx["<<i<<"]="<< meshx[i]<<" ";
    }
    output <<endl;
    for(int i=0;i<My;i++){
        output <<"meshy["<<i<<"]="<< meshy[i]<<" ";
    }
    output <<endl;

    output << "px="<<px<< " py="<<py<< endl;

    output.close();

    /*
    cerr << "***Meshing Stability Check***"<<endl;
    cerr << "Diffusion Check."<<endl;
    cerr << "sqrt(2*D_nm*tau) < Xstep => ";
    cerr << sqrt(2*D_nm*tau)<< " < " << MinMesh << endl;
    cerr << "Velocity Check."<<endl;
    cerr << "Velocity*tau < Xstep => ";
    cerr << MaxVelocity*tau << " < " << MinMesh <<endl;
    */
}


void CDmodel::CD_AddReceptors2D(){

    ReceptorArray = new Receptor[px];

    for(int m=0;m<Mx;m++){

        double a= xpin[m];

        for(int i=xb[m];i<xb[m+1]+1;i++){
            ReceptorArray[i].coordX=a+(i-xb[m])/meshx[m];
            ReceptorArray[i].AB=0;
            ReceptorArray[i].B=1e20;

        }
    }
}

void CDmodel::CD_Initialize(){

    #pragma omp parallel for
    for(int i=0;i<L;i++){
        CDmaterial[i].C_new=0;
        CDmaterial[i].C_old=0;
        CDmaterial[i].flag=0;
    }
}

void CDmodel::CD_InitialGuess2D(){

    CDmaterial=new CDElectrolyte [L];
    CD_Initialize();

    #pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            // concentration unit is M.
            int pointer = (px)*(j) + (i);

            // this can be 3D concentration (volume concentration)
            // the value 1 can has any unit. uM nM pM fM ...

            CDmaterial[pointer].C_new=C0;
            CDmaterial[pointer].C_old=C0;
        }
    }
}

void CDmodel::CD_Solver2D(){

    //using old calculate new, and then compare R, finally update to old.

    for(int j=1;j<TN+1;j++){

        #pragma omp parallel for
        for(int i=0;i<px;i++){
            for(int j=0;j<py;j++){

                //phi is new concentration
                //calculate concentration in next time step
                CD_inner2D(i,j);
            }
        }

        #pragma omp parallel for
        for(int i=0;i<px;i++){
            for(int j=0;j<py;j++){

                int pointer = (px)*(j) + (i);

                //copy new to old
                CDmaterial[pointer].C_old = CDmaterial[pointer].C_new;
            }
        }

        if(j%frame_step==0){

            stringstream name1, name2;
            string name11, name22;

            cerr <<"T="<<j*tau<<endl;
            name1<<"Concentration T="<<j*tau<<".txt";
            name11=name1.str();
            CD_PrintMaterial2D(name11.c_str());

            name2<<"Receptor T="<<j*tau<<".txt";
            name22=name2.str();
            CD_PrintReceptor2D(name22.c_str());

        }
        /*
        */
    }
}

void CDmodel::CD_inner2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);

    //double Rp=0, R=0;

    double x1=1,x2=1,x3=1,x4=1;
    double y1=1,y2=1,y3=1,y4=1;

    if (i==0){
        x2=0;
        x4=0;
        deltax=xstep_p;
        xstep_n=xstep_p;
        pointer_in=pointer;
    }
    if (i==px-1){
        x1=0;
        x3=0;
        deltax=xstep_n;
        xstep_p=xstep_n;
        pointer_ip=pointer;
    }
    if (j==0){
        y2=0;
        y4=0;
        deltay=ystep_p;
        ystep_n=ystep_p;
        pointer_jn=pointer;

        //Rp = -k1[A][B]+k-1[AB]
        //R=Rp*Area/Volume
        //Rp=(-1)*k_forward*CDmaterial[pointer].C_old*ReceptorArray[i].B+k_backward*ReceptorArray[i].AB; //surface density per second
        //R=(Rp)*deltax/deltax*deltay; //volume density per second
    }
    if (j==py-1){
        y1=0;
        y3=0;
        deltay=ystep_n;
        ystep_p=ystep_n;
        pointer_jp=pointer;
    }

    double Area=deltax*deltay;

    //F's unit is mole per second (because using FVM)
    double F =x1*(-1)*CFDmaterial[pointer].ui*Area*(CDmaterial[pointer_ip].C_old-CDmaterial[pointer].C_old)/xstep_p+x2*CFDmaterial[pointer].ui*Area*(CDmaterial[pointer_in].C_old-CDmaterial[pointer].C_old)/xstep_n
             +x3*D_nm*(CDmaterial[pointer_ip].C_old-CDmaterial[pointer].C_old)*deltay/xstep_p-x4*D_nm*(CDmaterial[pointer].C_old-CDmaterial[pointer_in].C_old)*deltay/xstep_n
             +y1*(-1)*CFDmaterial[pointer].vi*Area*(CDmaterial[pointer_jp].C_old-CDmaterial[pointer].C_old)/ystep_p+y2*CFDmaterial[pointer].vi*Area*(CDmaterial[pointer_jn].C_old-CDmaterial[pointer].C_old)/ystep_n
             +y3*D_nm*(CDmaterial[pointer_jp].C_old-CDmaterial[pointer].C_old)*deltax/ystep_p-y4*D_nm*(CDmaterial[pointer].C_old-CDmaterial[pointer_jn].C_old)*deltax/ystep_n;

    //A is volume density
    double A = CDmaterial[pointer].C_old+tau*(F)/Area;

    CDmaterial[pointer].C_new=A;

    if(A<0){
        cout << "Time resolution insufficient."<<endl;
        exit(0);
    }

    if(j==0){
        CDmaterial[pointer].C_new=0;
        ReceptorArray[i].AB=ReceptorArray[i].AB+A*Area/deltax;
    }else{
        CDmaterial[pointer].C_new=A;
    }

    // Binding Kinetics
    //Incomplete

    /*
    if(j==0){
        if(Rp<0){
            if(A*Area/deltax > ReceptorArray[i].B){
                if(abs(tau*Rp) <= ReceptorArray[i].B){
                    ReceptorArray[i].B=ReceptorArray[i].B+tau*Rp;
                    ReceptorArray[i].AB=ReceptorArray[i].AB-tau*Rp;
                    CDmaterial[pointer].C_new=A+tau*R;
                }else{
                    ReceptorArray[i].AB=ReceptorArray[i].AB+ReceptorArray[i].B;
                    CDmaterial[pointer].C_new=A-ReceptorArray[i].B*deltax/Area;
                    ReceptorArray[i].B=0;
                }
            }else{
                if(abs(tau*R) <= A){
                    ReceptorArray[i].B=ReceptorArray[i].B+tau*Rp;
                    ReceptorArray[i].AB=ReceptorArray[i].AB-tau*Rp;
                    CDmaterial[pointer].C_new=A+tau*R;
                }else{
                    ReceptorArray[i].B=ReceptorArray[i].B-A*Area/deltax;
                    ReceptorArray[i].AB=ReceptorArray[i].AB+A*Area/deltax;
                    CDmaterial[pointer].C_new=0;
                }
            }
        }else{
            if(Rp*tau<=ReceptorArray[i].AB){
                ReceptorArray[i].B=ReceptorArray[i].B+Rp*tau;
                ReceptorArray[i].AB=ReceptorArray[i].AB-Rp*tau;
                CDmaterial[pointer].C_new=A+Rp*tau*deltax/Area;
            }else{
                ReceptorArray[i].B=ReceptorArray[i].B+ReceptorArray[i].AB;
                CDmaterial[pointer].C_new=CDmaterial[pointer].C_new+ReceptorArray[i].AB*deltax/Area;
                ReceptorArray[i].AB=0;
            }
        }
    }
    */

}

void CDmodel::CD_PrintMaterial2D(string path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tY(2)\tC_new(3)\tC_old(4)\tui(5)\tvi(6)\tflag(7)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            output << mesh[pointer].coordX << '\t' << mesh[pointer].coordY << '\t'
                   <<CDmaterial[pointer].C_new << '\t' <<CDmaterial[pointer].C_old << '\t'
                   << CFDmaterial[pointer].ui << '\t' << CFDmaterial[pointer].vi << '\t' << CDmaterial[pointer].flag <<endl;

        }
    }

    output.close();
}

void CDmodel::CD_PrintReceptor2D(string path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tY(2)\tB(3)\tAB(4)\tflag(5)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        output << ReceptorArray[i].coordX << '\t' << ReceptorArray[i].coordY << '\t'
               << ReceptorArray[i].B << '\t' <<ReceptorArray[i].AB << '\t' <<ReceptorArray[i].flag<<endl;
    }

    output.close();
}
