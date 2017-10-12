#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#include "cdmodel.h"
#include "Parameter.h"


CDmodel::CDmodel()
{

}

CDmodel::~CDmodel()
{

}

void CDmodel::CDParameter(){

    t=50;       //time, s.
    tau=1e-3;  //time step0
    TN=t/tau+1;        //nubmer of steps

    AnalyteRadius_m=5e-9;         //radius [m]
    AnalyteRadius_nm=AnalyteRadius_m*1e9;

    eta_m=eta*1e18; //1e-3 [N*s/m^2]
    D=kb*Tamb/(6*M_PI*eta_m*AnalyteRadius_m);
    D_nm=D*1e18;

    /*A + B -> AB
     *
     *r1=k1[A][B]
     *r2=k2[AB]
     *
     *r2=r_-1
     *
     *R=(-k1[A][B]+k2[AB])*bottom area/Volume
     */

    k_forward=0;
    k_backward=0;

    k_forward1=1e3;
    k_backward1=1e-3;

    k_forward2=0;
    k_backward2=1e4;

}

void CDmodel::CDParameterSet2D(){

    CDParameter();
    CDStructureParameter2D();

    //Stability check
    double MinMesh=1e10;
    double MaxVelocity=300*1000;

    for(int i=0;i<Mx;i++){
        if(MinMesh>1/meshx[i]){
            MinMesh=1/meshx[i];
        }
    }

    for(int i=0;i<My;i++){
        if(MinMesh>1/meshy[i]){
            MinMesh=1/meshy[i];
        }
    }

    cerr << "***Meshing Stability Check***"<<endl;
    cerr << "Diffusion Check."<<endl;
    cerr << "sqrt(2*D_nm*tau) < Xstep => ";
    cerr << sqrt(2*D_nm*tau)<< " < " << MinMesh << endl;
    cerr << "Velocity Check."<<endl;
    cerr << "Velocity*tau < Xstep => ";
    cerr << MaxVelocity*tau << " < " << MinMesh <<endl;


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


}

void CDmodel::CDStructureParameter2D(){

    //Dimensions number
    DimensionTag=2;
    // xy in nm.
    ElectrodeGap=20*1e3;
    ElectrodeWidth=50*1e3;

    lx=400*1e3;
    ly=200*1e3;
    // xyz points iitial value
    px=py=1;

    // xyz block numbers.
    Mx=3;
    My=2;
    // set xy pins
    xpin=new double [Mx+1];
    ypin=new double [My+1];

    /*
    for(int i=0;i<Mx+1;i++){
        xpin[i]=0+lx/Mx*i;
    }
    */
    xpin[0]=0;
    xpin[1]=lx/2-ElectrodeGap/2-ElectrodeWidth*1;
    xpin[2]=lx/2+ElectrodeGap/2+ElectrodeWidth*1;
    xpin[3]=lx;
    /*
    for(int i=0;i<My+1;i++){
        ypin[i]=0+ly/My*i;
    }
    */
    ypin[0]=0;
    ypin[1]=ly/4;
    ypin[2]=ly;
    // set xy mesh steps
    meshx=new double [Mx];
    meshy=new double [My];

    // mesh unit = 1/nm = 1/step
    /*
    for(int i=0;i<Mx;i++){
        meshx[i]=1e-3*(i+1);
    }
    */

    meshx[0]=5e-4;
    meshx[1]=1e-3;
    meshx[2]=5e-4;

    /*
    for(int i=0;i<My;i++){
        meshy[i]=1e-3*(i+1);
    }
    */
    meshy[0]=1e-3;
    meshy[1]=5e-4;

    // points calculation
    for(int i=0;i<Mx;i++){
        px=px+meshx[i]*(xpin[i+1]-xpin[i]);
    }
    for(int i=0;i<My;i++){
        py=py+meshy[i]*(ypin[i+1]-ypin[i]);
    }
    // set xyz  point numbers till each block
    xb=new int [Mx+1];
    yb=new int [My+1];
    for(int i=1;i<Mx+1;i++){
        xb[0]=0;
        xb[i]=xb[i-1]+(xpin[i]-xpin[i-1])*meshx[i-1];
    }
    for(int i=1;i<My+1;i++){
        yb[0]=0;
        yb[i]=yb[i-1]+(ypin[i]-ypin[i-1])*meshy[i-1];
    }
    L1=px*py;
}

void CDmodel::CDAddReceptors2D(){

    ReceptorArray = new Receptor[px];

    for(int m=0;m<Mx;m++){

        double a= xpin[m];

        for(int i=xb[m];i<xb[m+1]+1;i++){
            ReceptorArray[i].coordX=a+(i-xb[m])/meshx[m];
            ReceptorArray[i].AB=0;
            ReceptorArray[i].B=1e20;

            /*
            //localized assignment
            if(ReceptorArray[i].coordX > 3*lx/8 && ReceptorArray[i].coordX < 5*lx/8 ){
                ReceptorArray[i].B=1e10;
            }
            */

        }
    }

}

void CDmodel::CDInitialGuess2D(){

    #pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){


            // concentration unit is M.
            int pointer = (px)*(j) + (i);

            // this can be 3D concentration (volume concentration)
            // the value 1 can has any unit. uM nM pM fM ...
            /*
            sample1[pointer].dop=1;
            sample1[pointer].phi=1;
            */

            if(pow(sample1[pointer].coordX-lx/2,2)+pow(sample1[pointer].coordY-ly/2,2) <= pow(20000,2)){
                sample1[pointer].dop=1;
                sample1[pointer].phi=1;
            }
            /*
            */
        }
    }

}

void CDmodel::CDSolver2D(){

    //using old calculate new, and then compare R, finally update to old.

    for(int j=1;j<TN;j++){

        #pragma omp parallel for
        for(int i=0;i<px;i++){
            for(int j=0;j<py;j++){

                //phi is new concentration
                //calculate concentration in next time step
                CDinner2D(i,j);
            }
        }

        #pragma omp parallel for
        for(int i=0;i<px;i++){
            for(int j=0;j<py;j++){


                int pointer = (px)*(j) + (i);

                //dop is old concentration
                //copy new to old, that is phi to dop. prepare for the next iteration.
                sample1[pointer].dop = sample1[pointer].phi;

            }
        }

        if(j%1000==0){

            stringstream name1, name2;
            string name11, name22;

            cerr <<"T="<<j*tau<<endl;
            name1<<"Concentration T="<<j*tau<<".txt";
            name11=name1.str();
            PrintMaterial2D(name11.c_str());

            name2<<"Receptor T="<<j*tau<<".txt";
            name22=name2.str();
            PrintReceptor2D(name22.c_str());

        }
        /*
        */
    }
}

void CDmodel::CDinner2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double deltax=abs(sample1[pointer_ip].coordX-sample1[pointer_in].coordX)/2;
    double deltay=abs(sample1[pointer_jp].coordY-sample1[pointer_jn].coordY)/2;
    double xstep_p=abs(sample1[pointer_ip].coordX-sample1[pointer].coordX);
    double xstep_n=abs(sample1[pointer_in].coordX-sample1[pointer].coordX);
    double ystep_p=abs(sample1[pointer_jp].coordY-sample1[pointer].coordY);
    double ystep_n=abs(sample1[pointer_jn].coordY-sample1[pointer].coordY);

    double Rp=0, R=0;

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
        Rp=(-1)*k1*sample1[pointer].dop*ReceptorArray[i].B+k2*ReceptorArray[i].AB; //surface density per second
        R=(Rp)*deltax/deltax*deltay; //volume density per second
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
    double F =x1*(-1)*sample1[pointer].ui*Area*(sample1[pointer_ip].dop-sample1[pointer].dop)/xstep_p+x2*sample1[pointer].ui*Area*(sample1[pointer_in].dop-sample1[pointer].dop)/xstep_n
             +x3*D_nm*(sample1[pointer_ip].dop-sample1[pointer].dop)*deltay/xstep_p-x4*D_nm*(sample1[pointer].dop-sample1[pointer_in].dop)*deltay/xstep_n
             +y1*(-1)*sample1[pointer].vi*Area*(sample1[pointer_jp].dop-sample1[pointer].dop)/ystep_p+y2*sample1[pointer].vi*Area*(sample1[pointer_jn].dop-sample1[pointer].dop)/ystep_n
             +y3*D_nm*(sample1[pointer_jp].dop-sample1[pointer].dop)*deltax/ystep_p-y4*D_nm*(sample1[pointer].dop-sample1[pointer_jn].dop)*deltax/ystep_n;

    //A is volume density
    double A = sample1[pointer].dop+tau*(F)/Area;

    sample1[pointer].phi=A;

    if(j==0){
        sample1[pointer].phi=0;
        ReceptorArray[i].AB=ReceptorArray[i].AB+A*Area/deltax;
    }else{
        sample1[pointer].phi=A;
    }

    //******************
    // Binding Kinetics
    //******************
    //Incomplete

    /*
    if(j==0){
        if(Rp<0){
            if(abs(tau*Rp) <= ReceptorArray[i].B){
                if(abs(tau*R) <= A){
                    ReceptorArray[i].B=ReceptorArray[i].B+tau*Rp;
                    ReceptorArray[i].AB=ReceptorArray[i].AB-tau*Rp;
                    sample1[pointer].phi=A+tau*R;
                }else{
                    ReceptorArray[i].B=ReceptorArray[i].B-A*Area/deltax;
                    ReceptorArray[i].AB=ReceptorArray[i].AB+A*Area/deltax;
                    sample1[pointer].phi=0;
                }
            }else{
                if(abs(tau*R) <= A){
                    ReceptorArray[i].AB=ReceptorArray[i].AB+ReceptorArray[i].B;
                    sample1[pointer].phi=A-ReceptorArray[i].B*deltax/Area;
                    ReceptorArray[i].B=0;
                }else{
                    double M=min(A*Area/deltax,ReceptorArray[i].B);
                    ReceptorArray[i].B=ReceptorArray[i].B-M;
                    ReceptorArray[i].AB=ReceptorArray[i].AB+M;
                    sample1[pointer].phi=A-M*deltax/Area;
                }
            }
        }else{
            if(Rp*tau<=ReceptorArray[i].AB){
                ReceptorArray[i].B=ReceptorArray[i].B+Rp*tau;
                ReceptorArray[i].AB=ReceptorArray[i].AB-Rp*tau;
                sample1[pointer].phi=A+Rp*tau*deltax/Area;
            }else{
                ReceptorArray[i].B=ReceptorArray[i].B+ReceptorArray[i].AB;
                sample1[pointer].phi=sample1[pointer].phi+ReceptorArray[i].AB*deltax/Area;
                ReceptorArray[i].AB=0;
            }
        }
    }else{
        sample1[pointer].phi=A;
    }
    */

}

void CDmodel::CDParameterSet3D(){

    CDParameter();
    CDStructureParameter3D();

    //Stability check
    double MinMesh=1e10;
    double MaxVelocity=300*1000;

    for(int i=0;i<Mx;i++){
        if(MinMesh>1/meshx[i]){
            MinMesh=1/meshx[i];
        }
    }

    for(int i=0;i<My;i++){
        if(MinMesh>1/meshy[i]){
            MinMesh=1/meshy[i];
        }
    }

    cerr << "***Meshing Stability Check***"<<endl;
    cerr << "Diffusion Check."<<endl;
    cerr << "sqrt(2*D_nm*tau) < Xstep => ";
    cerr << sqrt(2*D_nm*tau)<< " < " << MinMesh << endl;
    cerr << "Velocity Check."<<endl;
    cerr << "Velocity*tau < Xstep => ";
    cerr << MaxVelocity*tau << " < " << MinMesh <<endl;


    fstream output;
    output.open("Convection Diffusion Parameter", fstream::out | fstream::trunc);
    output << "lx="<<lx<< " ly="<<ly<< " lz="<<lz<< endl;
    output << "CC0="<<CC0<<endl;
    output << "SimTolPoisson="<<SimTolPoisson<<endl;
    output << "t="<<t<< " tau="<<tau<<endl;
    output << "AnalyteRadius_m="<<AnalyteRadius_m<<" mass="<<mass<<endl;
    output << "AnalyteValence="<<AnalyteValence<<" DotN="<<DotN<<endl;
    output << "Diffusion Constant="<<D<<" (m2/s) "<<D_nm<<" (nm2/s) "<<endl;
    output << "k1="<<k1<<" (M-1s-1), k2="<<k2<<" (s-1), ";
    output << "Ka="<<Ka<<", Kd="<<Kd<<endl;

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


}

void CDmodel::CDStructureParameter3D(){

    //Dimensions number
    DimensionTag=2;
    // xy in nm.
    lx=200*1000;
    ly=200*1000;
    // xyz points iitial value
    px=py=1;
    //BC
    LbX=0, UbX=lx, LbY=0, UbY=ly;
    // xyz block numbers.
    Mx=1;
    My=1;
    // set xy pins
    xpin=new double [Mx+1];
    ypin=new double [My+1];
    for(int i=0;i<Mx+1;i++){
        xpin[i]=0+lx/Mx*i;
    }
    /*
    xpin[0]=0;
    xpin[1]=lx/2-NWR-(1000-fmod(NWR,1000));
    xpin[2]=lx/2+NWR+(1000-fmod(NWR,1000));;
    xpin[3]=lx;
    */
    for(int i=0;i<My+1;i++){
        ypin[i]=0+ly/My*i;
    }
    /*
    ypin[0]=0;
    ypin[1]=50;
    ypin[2]=1000;
    ypin[3]=200000;
    */
    // set xy mesh steps
    meshx=new double [Mx];
    meshy=new double [My];
    // mesh unit = 1/nm = 1/step
    for(int i=0;i<Mx;i++){
        meshx[i]=0.001*(i+1);
    }

    /*
    meshx[0]=1e-3;
    meshx[1]=1e-2;
    meshx[2]=1e-3;
    */
    for(int i=0;i<My;i++){
        meshy[i]=0.001*(i+1);
    }
    /*
    meshx[3]=1e-2;
    meshx[4]=1e-3;

    meshy[0]=1;
    meshy[1]=2e-2;
    meshy[2]=1e-3;
    */
    // points calculation
    for(int i=0;i<Mx;i++){
        px=px+meshx[i]*(xpin[i+1]-xpin[i]);
    }
    for(int i=0;i<My;i++){
        py=py+meshy[i]*(ypin[i+1]-ypin[i]);
    }
    // set xyz  point numbers till each block
    xb=new int [Mx+1];
    yb=new int [My+1];
    for(int i=1;i<Mx+1;i++){
        xb[0]=0;
        xb[i]=xb[i-1]+(xpin[i]-xpin[i-1])*meshx[i-1];
    }
    for(int i=1;i<My+1;i++){
        yb[0]=0;
        yb[i]=yb[i-1]+(ypin[i]-ypin[i-1])*meshy[i-1];
    }
    L1=px*py;
}

void CDmodel::CDAddReceptors3D(){

    ReceptorArray = new Receptor[px];

    for(int m=0;m<Mx;m++){

        double a= xpin[m];

        for(int i=xb[m];i<xb[m+1]+1;i++){
            ReceptorArray[i].coordX=a+(i-xb[m])/meshx[m];
            ReceptorArray[i].B=1e10;
            ReceptorArray[i].AB=0;
        }
    }

}

void CDmodel::CDInitialGuess3D(){

    #pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){


            // concentration unit is M.
            int pointer = (px)*(j) + (i);

            sample1[pointer].dop=1;
            sample1[pointer].phi=1;
            /*
            if(sample1[pointer].coordX>3*lx/8 && sample1[pointer].coordX<5*lx/8 ){
                if(sample1[pointer].coordY>3*ly/8 && sample1[pointer].coordY<5*ly/8 ){
                }
            }
            */
        }
    }

}

void CDmodel::CDSolver3D(){

    for(int j=1;j<TN;j++){

        #pragma omp parallel for
        for(int i=0;i<px;i++){
            for(int j=0;j<py;j++){

                //phi is new concentration
                //calculate concentration in next time step
                CDinner3D(i,j);
            }
        }

        #pragma omp parallel for
        for(int i=0;i<px;i++){
            for(int j=0;j<py;j++){


                int pointer = (px)*(j) + (i);

                //dop is new concentration
                //update new to old
                sample1[pointer].dop = sample1[pointer].phi;

            }
        }

        if(j%1000==0){

            stringstream name1, name2;
            string name11, name22;

            name1<<"Concentration T="<<j*tau<<".txt";
            name11=name1.str();
            PrintMaterial2D(name11.c_str());

            name2<<"Receptor T="<<j*tau<<".txt";
            name22=name2.str();
            PrintReceptor2D(name22.c_str());

        }
        /*
        */
    }
}

void CDmodel::CDinner3D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double deltax=abs(sample1[pointer_ip].coordX-sample1[pointer_in].coordX)/2;
    double deltay=abs(sample1[pointer_jp].coordY-sample1[pointer_jn].coordY)/2;
    double xstep_p=abs(sample1[pointer_ip].coordX-sample1[pointer].coordX);
    double xstep_n=abs(sample1[pointer_in].coordX-sample1[pointer].coordX);
    double ystep_p=abs(sample1[pointer_jp].coordY-sample1[pointer].coordY);
    double ystep_n=abs(sample1[pointer_jn].coordY-sample1[pointer].coordY);

    double Rp=0, R=0;

    double x1=1,x2=1,x3=1,x4=1;
    double y1=1,y2=1,y3=1,y4=1;

    if (i==0){
        x2=0;
        x4=0;
        deltax=xstep_p;
        xstep_n=xstep_p;
    }
    if (i==px-1){
        x1=0;
        x3=0;
        deltax=xstep_n;
        xstep_p=xstep_n;
    }
    if (j==0){
        y2=0;
        y4=0;
        deltay=ystep_p;
        ystep_n=ystep_p;

        Rp=(-1)*k1*sample1[pointer].dop*ReceptorArray[i].B+k2*ReceptorArray[i].AB; //surface density per second
        R=(Rp)*deltax/deltax*deltay; //volume density per second
    }
    if (j==py-1){
        y1=0;
        y3=0;
        deltay=ystep_n;
        ystep_p=ystep_n;
    }

    double Area=deltax*deltay;

    //F is mole per second (because using FVM)
    double F =x1*(-1)*sample1[pointer].ui*Area*(sample1[pointer_ip].dop-sample1[pointer].dop)/xstep_p+x2*sample1[pointer].ui*Area*(sample1[pointer_in].dop-sample1[pointer].dop)/xstep_n
             +x3*D_nm*(sample1[pointer_ip].dop-sample1[pointer].dop)*deltay/xstep_p-x4*D_nm*(sample1[pointer].dop-sample1[pointer_in].dop)*deltay/xstep_n
             +y1*(-1)*sample1[pointer].vi*Area*(sample1[pointer_jp].dop-sample1[pointer].dop)/ystep_p+y2*sample1[pointer].vi*Area*(sample1[pointer_jn].dop-sample1[pointer].dop)/ystep_n
             +y3*D_nm*(sample1[pointer_jp].dop-sample1[pointer].dop)*deltax/ystep_p-y4*D_nm*(sample1[pointer].dop-sample1[pointer_jn].dop)*deltax/ystep_n;

    //A is volume density
    double A = sample1[pointer].dop+tau*(F)/Area;

    //old > new > R > update to old

    if(j==0 && (sample1[pointer].coordX > 3*lx/8 && sample1[pointer].coordX < 5*lx/8 ) ){
        if(Rp<0){
            if(abs(tau*Rp) <= ReceptorArray[i].B){
                if(abs(tau*R) <= A){
                    ReceptorArray[i].B=ReceptorArray[i].B+tau*Rp;
                    ReceptorArray[i].AB=ReceptorArray[i].AB-tau*Rp;
                    sample1[pointer].phi=A+tau*R;
                }else{
                    ReceptorArray[i].B=ReceptorArray[i].B-A*Area/deltax;
                    ReceptorArray[i].AB=ReceptorArray[i].AB+A*Area/deltax;
                    sample1[pointer].phi=0;
                }
            }else{
                if(abs(tau*R) <= A){
                    ReceptorArray[i].AB=ReceptorArray[i].AB+ReceptorArray[i].B;
                    sample1[pointer].phi=A-ReceptorArray[i].B*deltax/Area;
                    ReceptorArray[i].B=0;
                }else{
                    double M=min(A*Area/deltax,ReceptorArray[i].B);
                    ReceptorArray[i].B=ReceptorArray[i].B-M;
                    ReceptorArray[i].AB=ReceptorArray[i].AB+M;
                    sample1[pointer].phi=A-M*deltax/Area;
                }
            }
        }else{
            if(Rp*tau<=ReceptorArray[i].AB){
                ReceptorArray[i].B=ReceptorArray[i].B+Rp*tau;
                ReceptorArray[i].AB=ReceptorArray[i].AB-Rp*tau;
                sample1[pointer].phi=A+Rp*tau*deltax/Area;
            }else{
                ReceptorArray[i].B=ReceptorArray[i].B+ReceptorArray[i].AB;
                sample1[pointer].phi=sample1[pointer].phi+ReceptorArray[i].AB*deltax/Area;
                ReceptorArray[i].AB=0;
            }
        }
    }else{
        sample1[pointer].phi=A;
    }
}
