#include <cmath>
#include <fstream>
#include <sstream>

#include "ddmodel.h"
#include "Parameter.h"

DDmodel::DDmodel(){
}

DDmodel::~DDmodel(){
}

// DD model parameter setting
void DDmodel::DD_ParameterSet(){

    DD_BernoulliX();

    SimTolEC=SimTolHC=SimTolPoisson=1e-8;

    //voltage
    volS=0;
    volDe=1;
    volDs=0.1;
    volDi=0.1;
    volD=volDi;
    volB=0;
    volGe=2;
    volGs=0.1;
    volGi=0; //concomp -0.2
    volG=volGi;
    wfG=0;

    //tolerance
    DD_loop = 0;
    maxIter = 500;

    //Material
    Na=1e17;  // 1/cm3
    Nd=1e17;  // 1/cm3
    NaPlus=8e19;  // 1/cm3
    NdPlus=8e19;  // 1/cm3

    SiNc=2*pow(2*M_PI*Si_me*m0*kb*Tamb/pow(h,2),1.5);
    SiNv=2*pow(2*M_PI*Si_mh*m0*kb*Tamb/pow(h,2),1.5);
    Si_ni_m=sqrt(SiNc*SiNv)*exp((-1)*Si_Eg/(2*VT));
    Si_ni_cm=Si_ni_m*1e-6;
    Si_ni_nm=Si_ni_m*1e-27;

    ZnONc=2*pow(2*M_PI*ZnO_me*m0*kb*Tamb/pow(h,2),1.5);
    ZnONv=2*pow(2*M_PI*ZnO_mh*m0*kb*Tamb/pow(h,2),1.5);
    ZnO_ni_m=sqrt(ZnONc*ZnONv)*exp((-1)*ZnO_Eg/(2*VT));
    ZnO_ni_cm=ZnO_ni_m*1e-6;
    ZnO_ni_nm=ZnO_ni_m*1e-27;

    ni_m=Si_ni_m;
    ni_cm=Si_ni_cm;
    ni_nm=Si_ni_nm;

    Nai=Na/ni_cm;
    Ndi=Nd/ni_cm;
    NaPlusi=NaPlus/ni_cm;
    NdPlusi=NdPlus/ni_cm;

    mun_semi=0.12*1e18;
    mup_semi=0.12*1e18;
    mun_sub=0.12*1e18;
    mup_sub=0.12*1e18;

    //ion concentration
    //6.022*e-4 1/nm3 = 1mM = 1 mol/m3 = 1e-3 mol/L = 1e-3 M
    //9nm debye = 1mM,
    //3nm debye = 10mM
    //PH=7 => H+=1e-7 M = 1e-4*1e-3 M = 6.022*e-8 1/nm3
    CC0=1;
    C0=CC0*6.022e-4;
    mun_electrolyte=0.12*1e18;
    mup_electrolyte=0.12*1e18;

    DotNumber=5;
    AnalyteRadius=5;
    AnalyteValence=1;
    AnalytePermittivity=5;
    ReceptorLength=5;

}

void DDmodel::DD_NewAndInitialize(){

    DDmaterial=new Semiconductor [L];
    DD_Initialize();
}

void DDmodel::DD_InitialGuess2D(){


    switch(StructureFlag){

    case 1:
        DD_InitialGuessPNJunction2D();
        break;
    case 2:
        DD_InitialGuessMOSFET2D();
        break;
    case 3:
        DD_InitialGuessISFET2D();
        break;
    default:
        cout << "Undifined Device Structure @ DD_InitialGuess2D." << endl;
        exit(0);
    }

}

void DDmodel::DD_InitialGuessPNJunction2D(){

#pragma omp parallel for
    for(int i=0;i<L;i++){
       DDmaterial[i].k=Si_permi;
       DDmaterial[i].Type=1;
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer = (px)*(j) + (i);

            //setup P+ Drain
            //P
            DDmaterial[pointer].dop=-Nai;
            DDmaterial[pointer].phi=(volD-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
            DDmaterial[pointer].u=exp((-1)*volD/VT);
            DDmaterial[pointer].v=exp(volD/VT);
            DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
            DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
            DDmaterial[pointer].tau=DD_tauPCal(0);
            DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);

            //setup N+ Source
            if(mesh[pointer].coordX<lx/2){
                //N
                DDmaterial[pointer].dop=Ndi;
                DDmaterial[pointer].phi=(volS+VT*log(0.5*Ndi+sqrt(pow(0.5*Ndi,2)+1)));
                DDmaterial[pointer].u=exp((-1)*volS/VT);
                DDmaterial[pointer].v=exp(volS/VT);
                DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                DDmaterial[pointer].tau=DD_tauPCal(0);
                DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
            }
        }
    }
}

void DDmodel::DD_InitialGuessMOSFET2D(){

#pragma omp parallel for
    for(int i=0;i<L;i++){
       DDmaterial[i].k=Si_permi;
       DDmaterial[i].Type=1;
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer = (px)*(j) + (i);

            //setup P
            //P
            DDmaterial[pointer].dop=-Nai;
            DDmaterial[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
            DDmaterial[pointer].u=exp((-1)*volB/VT);
            DDmaterial[pointer].v=exp(volB/VT);
            DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
            DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
            DDmaterial[pointer].tau=DD_tauPCal(0);
            DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);

            //setup N+ Source
            if(mesh[pointer].coordX<=JunctionLength){
                if(mesh[pointer].coordY<=JunctionDepth){
                    //N+
                    DDmaterial[pointer].dop=NdPlusi;
                    DDmaterial[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                    DDmaterial[pointer].u=exp((-1)*volS/VT);
                    DDmaterial[pointer].v=exp(volS/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                    DDmaterial[pointer].tau=DD_tauPCal(0);
                    DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                }
            }

            //setup N+ Drain
            if(mesh[pointer].coordX>=lx-JunctionLength){
                if(mesh[pointer].coordY<=JunctionDepth){
                    //N+
                    DDmaterial[pointer].dop=NdPlusi;
                    DDmaterial[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                    DDmaterial[pointer].u=exp((-1)*volD/VT);
                    DDmaterial[pointer].v=exp(volD/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                    DDmaterial[pointer].tau=DD_tauPCal(0);
                    DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                }
            }
        }
    }
}

void DDmodel::DD_InitialGuessISFET2D(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer = (px)*(j) + (i);

            //substrate
            if(mesh[pointer].coordY<=SubstrateThickness ){
                //P
                DDmaterial[pointer].Type=1;
                DDmaterial[pointer].k=Si_permi;
                DDmaterial[pointer].dop=-Nai;
                DDmaterial[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                DDmaterial[pointer].u=exp((-1)*volB/VT);
                DDmaterial[pointer].v=exp(volB/VT);
                DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                DDmaterial[pointer].tau=DD_tauPCal(0);
                DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
            }

            //SD
            if(mesh[pointer].coordY>SubstrateThickness-JunctionDepth && mesh[pointer].coordY<=SubstrateThickness ){

                if(mesh[pointer].coordX<JunctionLength){
                    //N+
                    DDmaterial[pointer].Type=1;
                    DDmaterial[pointer].k=Si_permi;
                    DDmaterial[pointer].dop=NdPlusi;
                    DDmaterial[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                    DDmaterial[pointer].u=exp((-1)*volS/VT);
                    DDmaterial[pointer].v=exp(volS/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                    DDmaterial[pointer].tau=DD_tauPCal(0);
                    DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);

                }
                if(mesh[pointer].coordX>lx-JunctionLength){
                    //N+
                    DDmaterial[pointer].Type=1;
                    DDmaterial[pointer].k=Si_permi;
                    DDmaterial[pointer].dop=NdPlusi;
                    DDmaterial[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                    DDmaterial[pointer].u=exp((-1)*volD/VT);
                    DDmaterial[pointer].v=exp(volD/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                    DDmaterial[pointer].tau=DD_tauPCal(0);
                    DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                }
            }

            //oxide
            if(mesh[pointer].coordY>SubstrateThickness && mesh[pointer].coordY<=SubstrateThickness+Tox ){
                DDmaterial[pointer].Type=2;
                DDmaterial[pointer].k=SiO2_permi;
                DDmaterial[pointer].phi=volG;
                DDmaterial[pointer].dop=0;
                DDmaterial[pointer].u=exp((-1)*volG/VT);
                DDmaterial[pointer].v=exp(volG/VT);
                DDmaterial[pointer].mun=0;
                DDmaterial[pointer].mup=0;
                DDmaterial[pointer].tau=0;
                DDmaterial[pointer].r=0;
            }

            //electrolyte
            if(mesh[pointer].coordY>SubstrateThickness+Tox){
                DDmaterial[pointer].Type=3;
                DDmaterial[pointer].k=Water_permi;
                DDmaterial[pointer].phi=volG;
                DDmaterial[pointer].dop=0;
                DDmaterial[pointer].u=exp((-1)*volG/VT);
                DDmaterial[pointer].v=exp(volG/VT);
                DDmaterial[pointer].mun=0.12*1e18;
                DDmaterial[pointer].mup=0.12*1e18;
                DDmaterial[pointer].tau=DD_tauNCal(Si_ni_cm);
                DDmaterial[pointer].r=0;
            }

        }
    }

    //find boundary for sub:tox, tox:eletrolyte
    for(int j=0;j<py-1;j++){

        int pointer = (px)*(j) + (px/2);
        int pointer_jn = (px)*(j-1) + (px/2);
        int pointer_jp = (px)*(j+1) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_jp].Type==2){
            SubstrateTOP=j;
        }

        if(DDmaterial[pointer].Type==3 && DDmaterial[pointer_jn].Type==2){
            ElectrolyteBottom=j;
        }
    }

    if(SubstrateTOP==0||ElectrolyteBottom==0){
        cout << "Boundary undifined."<<endl;
        cout << "SubstrateTOP=" <<SubstrateTOP <<endl;
        cout << "ElectrolyteBottom="<<ElectrolyteBottom<<endl;
        exit(0);
    }
}

void DDmodel::DD_InitialGuess3D(){

    DDmaterial=new Semiconductor [L];
    DD_Initialize();

    switch(StructureFlag){

    case 1:
        DD_InitialGuessPNJunction3D();
        break;
    case 2:
        DD_InitialGuessMOSFET3D();
        break;
    case 3:
        DD_InitialGuessISFET3D();
        break;
    case 4:
        DD_InitialGuess1NWR3D();
        DD_FindNWRBC1();
        break;
    case 5:
        DD_InitialGuess2NWR3D();
        DD_FindNWRBC1();
        DD_FindNWRBC2();
        break;
    default:
        cout << "Undifined Device Structure @ DD_InitialGuess3D." << endl;
        exit(0);
    }

}

void DDmodel::DD_InitialGuessPNJunction3D(){

    StructureFlag=1;

#pragma omp parallel for
    for(int i=0;i<L;i++){
       DDmaterial[i].k=Si_permi;
       DDmaterial[i].Type=1;
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //setup P Drain
                //P
                DDmaterial[pointer].dop=-Nai;
                DDmaterial[pointer].phi=(volD-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                DDmaterial[pointer].u=exp((-1)*volD/VT);
                DDmaterial[pointer].v=exp(volD/VT);
                DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                DDmaterial[pointer].tau=DD_tauPCal(0);
                DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);

                //setup N Source
                if(mesh[pointer].coordX<lx/2){
                    //N
                    DDmaterial[pointer].dop=Ndi;
                    DDmaterial[pointer].phi=(volS+VT*log(0.5*Ndi+sqrt(pow(0.5*Ndi,2)+1)));
                    DDmaterial[pointer].u=exp((-1)*volS/VT);
                    DDmaterial[pointer].v=exp(volS/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                    DDmaterial[pointer].tau=DD_tauPCal(0);
                    DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                }
            }
        }
    }
}

void DDmodel::DD_InitialGuessMOSFET3D(){

    //StructureFlag=2;

#pragma omp parallel for
    for(int i=0;i<L;i++){
       DDmaterial[i].k=Si_permi;
       DDmaterial[i].Type=1;
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //setup P
                //P
                DDmaterial[pointer].dop=-Nai;
                DDmaterial[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                DDmaterial[pointer].u=exp((-1)*volB/VT);
                DDmaterial[pointer].v=exp(volB/VT);
                DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                DDmaterial[pointer].tau=DD_tauPCal(0);
                DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);

                //setup N+ Source
                if(mesh[pointer].coordX<=JunctionLength){
                    if(mesh[pointer].coordZ<=JunctionDepth){
                        //N+
                        DDmaterial[pointer].dop=NdPlusi;
                        DDmaterial[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                        DDmaterial[pointer].u=exp((-1)*volS/VT);
                        DDmaterial[pointer].v=exp(volS/VT);
                        DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                        DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                        DDmaterial[pointer].tau=DD_tauPCal(0);
                        DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                    }
                }

                //setup N+ Drain
                if(mesh[pointer].coordX>=lx-JunctionLength){
                    if(mesh[pointer].coordZ<=JunctionDepth){
                        //N+
                        DDmaterial[pointer].dop=NdPlusi;
                        DDmaterial[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                        DDmaterial[pointer].u=exp((-1)*volD/VT);
                        DDmaterial[pointer].v=exp(volD/VT);
                        DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                        DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                        DDmaterial[pointer].tau=DD_tauPCal(0);
                        DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                    }
                }
            }
        }
    }
}

void DDmodel::DD_InitialGuessISFET3D(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //substrate
                if(mesh[pointer].coordZ<=SubstrateThickness ){
                    //P
                    DDmaterial[pointer].Type=1;
                    DDmaterial[pointer].k=Si_permi;
                    DDmaterial[pointer].dop=-Nai;
                    DDmaterial[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                    DDmaterial[pointer].u=exp((-1)*volB/VT);
                    DDmaterial[pointer].v=exp(volB/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                    DDmaterial[pointer].tau=DD_tauPCal(0);
                    DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                }

                //SD
                if(mesh[pointer].coordZ>SubstrateThickness-JunctionDepth && mesh[pointer].coordZ<=SubstrateThickness ){

                    if(mesh[pointer].coordX<JunctionLength){
                        //N+
                        DDmaterial[pointer].Type=1;
                        DDmaterial[pointer].k=Si_permi;
                        DDmaterial[pointer].dop=NdPlusi;
                        DDmaterial[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                        DDmaterial[pointer].u=exp((-1)*volS/VT);
                        DDmaterial[pointer].v=exp(volS/VT);
                        DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                        DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                        DDmaterial[pointer].tau=DD_tauPCal(0);
                        DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                    }
                    if(mesh[pointer].coordX>lx-JunctionLength){
                        //N+
                        DDmaterial[pointer].Type=1;
                        DDmaterial[pointer].k=Si_permi;
                        DDmaterial[pointer].dop=NdPlusi;
                        DDmaterial[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                        DDmaterial[pointer].u=exp((-1)*volD/VT);
                        DDmaterial[pointer].v=exp(volD/VT);
                        DDmaterial[pointer].mun=DD_munCal(Tamb, 0, 1); // max Na Nd
                        DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, 1);
                        DDmaterial[pointer].tau=DD_tauPCal(0);
                        DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);
                    }
                }

                //oxide
                if(mesh[pointer].coordZ>SubstrateThickness && mesh[pointer].coordZ<=SubstrateThickness+Tox ){
                    DDmaterial[pointer].Type=2;
                    DDmaterial[pointer].k=SiO2_permi;
                    DDmaterial[pointer].phi=volG;
                    DDmaterial[pointer].dop=0;
                    DDmaterial[pointer].u=exp((-1)*volG/VT);
                    DDmaterial[pointer].v=exp(volG/VT);
                    DDmaterial[pointer].mun=0;
                    DDmaterial[pointer].mup=0;
                    DDmaterial[pointer].tau=0;
                    DDmaterial[pointer].r=0;
                }

                //electrolyte
                if(mesh[pointer].coordZ>SubstrateThickness+Tox){
                    DDmaterial[pointer].Type=3;
                    DDmaterial[pointer].k=Water_permi;
                    DDmaterial[pointer].phi=volG;
                    DDmaterial[pointer].dop=0;
                    DDmaterial[pointer].u=exp((-1)*volG/VT);
                    DDmaterial[pointer].v=exp(volG/VT);
                    DDmaterial[pointer].mun=0.12*1e18;
                    DDmaterial[pointer].mup=0.12*1e18;
                    DDmaterial[pointer].tau=DD_tauNCal(Si_ni_cm);
                    DDmaterial[pointer].r=0;
                }
            }
        }
    }

    //find boundary for sub:tox, tox:eletrolyte
    for(int k=0;k<pz-1;k++){

        int pointer = (py)*(px)*(k) + (px)*(py/2) + (px/2);
        int pointer_kn = (py)*(px)*(k-1) + (px)*(py/2) + (px/2);
        int pointer_kp = (py)*(px)*(k+1) + (px)*(py/2) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_kp].Type==2){
            SubstrateTOP=k;
            //cerr << DDmaterial[pointer_jp].Type <<'\t'<<j<<endl;
        }

        if(DDmaterial[pointer].Type==3 && DDmaterial[pointer_kn].Type==2){
            ElectrolyteBottom=k;
            //cerr << DDmaterial[pointer_jn].Type <<'\t'<<j<<endl;
        }
    }

    if(SubstrateTOP==0||ElectrolyteBottom==0){
        cerr << "Boundary undifined."<<endl;
        cerr << "SubstrateTOP=" <<SubstrateTOP <<endl;
        cerr << "ElectrolyteBottom="<<ElectrolyteBottom<<endl;
    }
}

void DDmodel::DD_InitialGuess1NWR3D(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //substrate
                if(mesh[pointer].coordZ<=SubstrateThickness ){
                    DDmaterial[pointer].Type=4;
                    DDmaterial[pointer].k=Si_permi;
                    DDmaterial[pointer].dop=0;
                    DDmaterial[pointer].phi=volB;
                    DDmaterial[pointer].u=exp((-1)*volB/VT);
                    DDmaterial[pointer].v=exp(volB/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, DDmaterial[pointer].Type); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, DDmaterial[pointer].Type);
                    DDmaterial[pointer].tau=DD_tauPCal(0);
                    DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                    DDmaterial[pointer].rho=0;
                }

                //electrolyte
                if(mesh[pointer].coordZ>SubstrateThickness+BOX){
                    DDmaterial[pointer].Type=3;
                    DDmaterial[pointer].k=Water_permi;
                    //n+
                    //DDmaterial[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2-VT*log(0.5*Ndi_SD+sqrt(pow(0.5*Ndi_SD,2)+1)))+volG;
                    //P+
                    //DDmaterial[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2+VT*log(0.5*Nai_SD+sqrt(pow(0.5*Nai_SD,2)+1)))+volG;
                    //metal
                    //phi_new[pointer]=(Si_chi+(Si_Eg)/2)-(wfG)+volG;
                    //Intrinsic
                    DDmaterial[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2)+volG;
                    DDmaterial[pointer].dop=0;
                    DDmaterial[pointer].u=exp((-1)*volG/VT);
                    DDmaterial[pointer].v=exp(volG/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, DDmaterial[pointer].Type); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, DDmaterial[pointer].Type);
                    DDmaterial[pointer].tau=DD_tauNCal(Si_ni_cm);
                    DDmaterial[pointer].r=0;
                    DDmaterial[pointer].rho=0;
                }

                //box oxide
                //Protruding
                //if(mesh[pointer].coordZ>SubstrateThickness && mesh[pointer].coordZ<=SubstrateThickness+BOX ){
                //Burried
                if(mesh[pointer].coordZ>SubstrateThickness && mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz+Tox ){
                    DDmaterial[pointer].Type=2;
                    DDmaterial[pointer].k=SiO2_permi;
                    DDmaterial[pointer].phi=volB;
                    DDmaterial[pointer].dop=0;
                    DDmaterial[pointer].u=0;
                    DDmaterial[pointer].v=0;
                    DDmaterial[pointer].mun=0;
                    DDmaterial[pointer].mup=0;
                    DDmaterial[pointer].tau=0;
                    DDmaterial[pointer].r=0;
                    DDmaterial[pointer].rho=0;
                }

                //NWR oxide
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz+Tox && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery1+NWRradiusy+Tox && mesh[pointer].coordY>=lx/2-NWRradiusy-Tox){
                        DDmaterial[pointer].Type=2;
                        DDmaterial[pointer].k=SiO2_permi;
                        DDmaterial[pointer].phi=volB;
                        DDmaterial[pointer].dop=0;
                        DDmaterial[pointer].u=0;
                        DDmaterial[pointer].v=0;
                        DDmaterial[pointer].mun=0;
                        DDmaterial[pointer].mup=0;
                        DDmaterial[pointer].tau=0;
                        DDmaterial[pointer].r=0;
                        DDmaterial[pointer].rho=0;
                    }
                }

                //NWR1 channel
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery1+NWRradiusy && mesh[pointer].coordY>=NWRCentery1-NWRradiusy){
                        //P+
                        DDmaterial[pointer].Type=1;
                        DDmaterial[pointer].k=Si_permi;
                        DDmaterial[pointer].dop=-Nai;
                        DDmaterial[pointer].phi=((volD+volS)-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                        DDmaterial[pointer].u=exp((-1)*(volD+volS)/VT);
                        DDmaterial[pointer].v=exp((volD+volS)/VT);
                        DDmaterial[pointer].mun=DD_munCal(Tamb, Na, DDmaterial[pointer].Type); // max Na Nd
                        DDmaterial[pointer].mup=DD_mupCal(Tamb, Na, DDmaterial[pointer].Type);
                        DDmaterial[pointer].tau=DD_tauPCal(0);
                        DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                        DDmaterial[pointer].rho=0;

                        if(mesh[pointer].coordX<JunctionLength){
                            DDmaterial[pointer].Type=1;
                            DDmaterial[pointer].k=Si_permi;
                            DDmaterial[pointer].dop=-NaPlusi;
                            DDmaterial[pointer].phi=(volS-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            DDmaterial[pointer].u=exp((-1)*volS/VT);
                            DDmaterial[pointer].v=exp(volS/VT);
                            DDmaterial[pointer].mun=DD_munCal(Tamb, NaPlus, DDmaterial[pointer].Type); // max Na Nd
                            DDmaterial[pointer].mup=DD_mupCal(Tamb, NaPlus, DDmaterial[pointer].Type);
                            DDmaterial[pointer].tau=DD_tauPCal(0);
                            DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                            DDmaterial[pointer].rho=0;

                        }
                        if(mesh[pointer].coordX>ly-JunctionLength){
                            DDmaterial[pointer].Type=1;
                            DDmaterial[pointer].k=Si_permi;
                            DDmaterial[pointer].dop=-NaPlusi;
                            DDmaterial[pointer].phi=(volD-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            DDmaterial[pointer].u=exp((-1)*volD/VT);
                            DDmaterial[pointer].v=exp(volD/VT);
                            DDmaterial[pointer].mun=DD_munCal(Tamb, NaPlus, DDmaterial[pointer].Type); // max Na Nd
                            DDmaterial[pointer].mup=DD_mupCal(Tamb, NaPlus, DDmaterial[pointer].Type);
                            DDmaterial[pointer].tau=DD_tauPCal(0);
                            DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                            DDmaterial[pointer].rho=0;
                        }
                    }
                }
            }
        }
    }

}

void DDmodel::DD_InitialGuess2NWR3D(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //substrate
                if(mesh[pointer].coordZ<=SubstrateThickness ){
                    DDmaterial[pointer].Type=4;
                    DDmaterial[pointer].k=Si_permi;
                    DDmaterial[pointer].dop=0;
                    DDmaterial[pointer].phi=volB;
                    DDmaterial[pointer].u=exp((-1)*volB/VT);
                    DDmaterial[pointer].v=exp(volB/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, DDmaterial[pointer].Type); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, DDmaterial[pointer].Type);
                    DDmaterial[pointer].tau=DD_tauPCal(0);
                    DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                    DDmaterial[pointer].rho=0;
                }

                //electrolyte
                if(mesh[pointer].coordZ>SubstrateThickness+BOX){
                    DDmaterial[pointer].Type=3;
                    DDmaterial[pointer].k=Water_permi;
                    //n+
                    //DDmaterial[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2-VT*log(0.5*Ndi_SD+sqrt(pow(0.5*Ndi_SD,2)+1)))+volG;
                    //P+
                    //DDmaterial[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2+VT*log(0.5*Nai_SD+sqrt(pow(0.5*Nai_SD,2)+1)))+volG;
                    //metal
                    //phi_new[pointer]=(Si_chi+(Si_Eg)/2)-(wfG)+volG;
                    //Intrinsic
                    DDmaterial[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2)+volG;
                    DDmaterial[pointer].dop=0;
                    DDmaterial[pointer].u=exp((-1)*volG/VT);
                    DDmaterial[pointer].v=exp(volG/VT);
                    DDmaterial[pointer].mun=DD_munCal(Tamb, 0, DDmaterial[pointer].Type); // max Na Nd
                    DDmaterial[pointer].mup=DD_mupCal(Tamb, 0, DDmaterial[pointer].Type);
                    DDmaterial[pointer].tau=DD_tauNCal(Si_ni_cm);
                    DDmaterial[pointer].r=0;
                    DDmaterial[pointer].rho=0;
                }

                //box oxide
                //Burried
                if(mesh[pointer].coordZ>SubstrateThickness && mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz+Tox ){
                    DDmaterial[pointer].Type=2;
                    DDmaterial[pointer].k=SiO2_permi;
                    DDmaterial[pointer].phi=volB;
                    DDmaterial[pointer].dop=0;
                    DDmaterial[pointer].u=0;
                    DDmaterial[pointer].v=0;
                    DDmaterial[pointer].mun=0;
                    DDmaterial[pointer].mup=0;
                    DDmaterial[pointer].tau=0;
                    DDmaterial[pointer].r=0;
                    DDmaterial[pointer].rho=0;
                }

                //NWR oxide
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz+Tox && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery1+NWRradiusy+Tox && mesh[pointer].coordY>=lx/2-NWRradiusy-Tox){
                        DDmaterial[pointer].Type=2;
                        DDmaterial[pointer].k=SiO2_permi;
                        DDmaterial[pointer].phi=volB;
                        DDmaterial[pointer].dop=0;
                        DDmaterial[pointer].u=0;
                        DDmaterial[pointer].v=0;
                        DDmaterial[pointer].mun=0;
                        DDmaterial[pointer].mup=0;
                        DDmaterial[pointer].tau=0;
                        DDmaterial[pointer].r=0;
                        DDmaterial[pointer].rho=0;
                    }
                }

                //NWR1 channel
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery1+NWRradiusy && mesh[pointer].coordY>=NWRCentery1-NWRradiusy){
                        //P+
                        DDmaterial[pointer].Type=1;
                        DDmaterial[pointer].k=Si_permi;
                        DDmaterial[pointer].dop=-Nai;
                        DDmaterial[pointer].phi=((volD+volS)-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                        DDmaterial[pointer].u=exp((-1)*(volD+volS)/VT);
                        DDmaterial[pointer].v=exp((volD+volS)/VT);
                        DDmaterial[pointer].mun=DD_munCal(Tamb, Na, DDmaterial[pointer].Type); // max Na Nd
                        DDmaterial[pointer].mup=DD_mupCal(Tamb, Na, DDmaterial[pointer].Type);
                        DDmaterial[pointer].tau=DD_tauPCal(0);
                        DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                        DDmaterial[pointer].rho=0;

                        if(mesh[pointer].coordX<JunctionLength){
                            DDmaterial[pointer].Type=1;
                            DDmaterial[pointer].k=Si_permi;
                            DDmaterial[pointer].dop=-NaPlusi;
                            DDmaterial[pointer].phi=(volS-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            DDmaterial[pointer].u=exp((-1)*volS/VT);
                            DDmaterial[pointer].v=exp(volS/VT);
                            DDmaterial[pointer].mun=DD_munCal(Tamb, NaPlus, DDmaterial[pointer].Type); // max Na Nd
                            DDmaterial[pointer].mup=DD_mupCal(Tamb, NaPlus, DDmaterial[pointer].Type);
                            DDmaterial[pointer].tau=DD_tauPCal(0);
                            DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                            DDmaterial[pointer].rho=0;

                        }
                        if(mesh[pointer].coordX>ly-JunctionLength){
                            DDmaterial[pointer].Type=1;
                            DDmaterial[pointer].k=Si_permi;
                            DDmaterial[pointer].dop=-NaPlusi;
                            DDmaterial[pointer].phi=(volD-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            DDmaterial[pointer].u=exp((-1)*volD/VT);
                            DDmaterial[pointer].v=exp(volD/VT);
                            DDmaterial[pointer].mun=DD_munCal(Tamb, NaPlus, DDmaterial[pointer].Type); // max Na Nd
                            DDmaterial[pointer].mup=DD_mupCal(Tamb, NaPlus, DDmaterial[pointer].Type);
                            DDmaterial[pointer].tau=DD_tauPCal(0);
                            DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                            DDmaterial[pointer].rho=0;
                        }
                    }
                }

                //NWR2 channel
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery2+NWRradiusy && mesh[pointer].coordY>=NWRCentery2-NWRradiusy){
                        //P+
                        DDmaterial[pointer].Type=1;
                        DDmaterial[pointer].k=Si_permi;
                        DDmaterial[pointer].dop=-Nai;
                        DDmaterial[pointer].phi=((volD+volS)-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                        DDmaterial[pointer].u=exp((-1)*(volD+volS)/VT);
                        DDmaterial[pointer].v=exp((volD+volS)/VT);
                        DDmaterial[pointer].mun=DD_munCal(Tamb, Na, DDmaterial[pointer].Type); // max Na Nd
                        DDmaterial[pointer].mup=DD_mupCal(Tamb, Na, DDmaterial[pointer].Type);
                        DDmaterial[pointer].tau=DD_tauPCal(0);
                        DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                        DDmaterial[pointer].rho=0;

                        if(mesh[pointer].coordX<JunctionLength){
                            DDmaterial[pointer].Type=1;
                            DDmaterial[pointer].k=Si_permi;
                            DDmaterial[pointer].dop=-NaPlusi;
                            DDmaterial[pointer].phi=(volS-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            DDmaterial[pointer].u=exp((-1)*volS/VT);
                            DDmaterial[pointer].v=exp(volS/VT);
                            DDmaterial[pointer].mun=DD_munCal(Tamb, NaPlus, DDmaterial[pointer].Type); // max Na Nd
                            DDmaterial[pointer].mup=DD_mupCal(Tamb, NaPlus, DDmaterial[pointer].Type);
                            DDmaterial[pointer].tau=DD_tauPCal(0);
                            DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                            DDmaterial[pointer].rho=0;

                        }
                        if(mesh[pointer].coordX>ly-JunctionLength){
                            DDmaterial[pointer].Type=1;
                            DDmaterial[pointer].k=Si_permi;
                            DDmaterial[pointer].dop=-NaPlusi;
                            DDmaterial[pointer].phi=(volD-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            DDmaterial[pointer].u=exp((-1)*volD/VT);
                            DDmaterial[pointer].v=exp(volD/VT);
                            DDmaterial[pointer].mun=DD_munCal(Tamb, NaPlus, DDmaterial[pointer].Type); // max Na Nd
                            DDmaterial[pointer].mup=DD_mupCal(Tamb, NaPlus, DDmaterial[pointer].Type);
                            DDmaterial[pointer].tau=DD_tauPCal(0);
                            DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);
                            DDmaterial[pointer].rho=0;
                        }
                    }
                }
            }
        }
    }

}

double DDmodel::DD_PoissonSolver2D(){

    DD_loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        DD_loop++;

        errPhi=DD_PoissonGaussSeidel2D();

        if(errPhi_max < errPhi) {errPhi_max=errPhi;}

        if(DD_loop%1000==0)
        cout <<"PS:"<< DD_loop <<"\t" <<errPhi<<"\t"<<errPhi_max<<endl;

        if(DD_loop%100000==0)
        DD_PrintMaterial2D("Poisson_temp.txt");

    }while(errPhi>SimTolPoisson);

    return errPhi_max;

}

double DDmodel::DD_PoissonGaussSeidel2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double phik=DDmaterial[pointer].phi;

            DDmaterial[pointer].phi=DD_PoissonGaussSeidelInner2D(i,j);

            DDmaterial[pointer].r=DD_SRHrecomb2D(i,j);

            double error=abs(DDmaterial[pointer].phi-phik);

            error=error/(abs(phik)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    DD_PoissonBC2D();

    return max_val;

}

double DDmodel::DD_PoissonGaussSeidelInner2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double permitivity_ip=DDmaterial[pointer].k*DDmaterial[pointer_ip].k / (0.5*DDmaterial[pointer_ip].k+0.5*DDmaterial[pointer].k);
    double permitivity_in=DDmaterial[pointer].k*DDmaterial[pointer_in].k / (0.5*DDmaterial[pointer_in].k+0.5*DDmaterial[pointer].k);
    double permitivity_jp=DDmaterial[pointer].k*DDmaterial[pointer_jp].k / (0.5*DDmaterial[pointer_jp].k+0.5*DDmaterial[pointer].k);
    double permitivity_jn=DDmaterial[pointer].k*DDmaterial[pointer_jn].k / (0.5*DDmaterial[pointer_jn].k+0.5*DDmaterial[pointer].k);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);

    double f,df,phik;
    double volume=deltax*deltay;

    phik=DDmaterial[pointer].phi;

    int Temp=DDmaterial[pointer].Type;

    switch(Temp){

    case 1: //1=channel
        f=volume*ni_nm*(DDmaterial[pointer].u*exp(phik/VT)-DDmaterial[pointer].v*exp((-1)*phik/VT)-DDmaterial[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(DDmaterial[pointer].u*exp(phik/VT)+DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 2: //2=insulator
        f=volume*DDmaterial[pointer].rho/e0;
        df=0;
        break;
    case 3: //3=electrolyte
        f=volume*C0*(DDmaterial[pointer].u*exp(phik/VT)-DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0;
       df=volume*C0*(DDmaterial[pointer].u*exp(phik/VT)+DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 4: //4=Substrate
        f=volume*ni_nm*(DDmaterial[pointer].u*exp(phik/VT)-DDmaterial[pointer].v*exp((-1)*phik/VT)-DDmaterial[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(DDmaterial[pointer].u*exp(phik/VT)+DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 5: //5=Analyte
        f=volume*DDmaterial[pointer].rho/e0;
        df=0;
        break;
    default:
        cout << "Error! Undefined Type for Poisson Solver @ PoissonGaussSeidelInner2D."<<endl;
        exit(0);
    }

    return (((permitivity_ip*DDmaterial[pointer_ip].phi/xstep_p+permitivity_in*DDmaterial[pointer_in].phi/xstep_n)*deltay
            +(permitivity_jp*DDmaterial[pointer_jp].phi/ystep_p+permitivity_jn*DDmaterial[pointer_jn].phi/ystep_n)*deltax
            + f - df*phik )
            /
            ((permitivity_ip/xstep_p+permitivity_in/xstep_n)*deltay
            +(permitivity_jp/ystep_p+permitivity_jn/ystep_n)*deltax - df));
}

double DDmodel::DD_PoissonSolver3D(){

    DD_loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        DD_loop++;

        errPhi=DD_PoissonGaussSeidel3D();

        if(errPhi_max < errPhi) {errPhi_max=errPhi;}

        if(DD_loop%1000==0)
        cerr <<"PS:"<< DD_loop <<"\t" <<errPhi<<"\t"<<errPhi_max<<endl;

        /*
        if(DD_loop%10000==0){
            RhoCalculation3D();
            PrintMaterial3D("Poisson_temp.txt");
        }
        */

    }while(errPhi>SimTolPoisson);

    return errPhi_max;

}

double DDmodel::DD_PoissonGaussSeidel3D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<pz-1; k++) {

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                double phik=DDmaterial[pointer].phi;

                DDmaterial[pointer].phi=DD_PoissonGaussSeidelInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=abs(DDmaterial[pointer].phi-phik);

                error=error/(abs(phik)+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_PoissonBC3D();

    return max_val;

}

double DDmodel::DD_PoissonGaussSeidelInner3D(int i, int j, int k){

    int pointer = (py)*(px)*(k) + (px)*(j) + (i);
    int pointer_ip = (py)*(px)*(k) + (px)*(j) + (i+1);
    int pointer_in = (py)*(px)*(k) + (px)*(j) + (i-1);
    int pointer_jp = (py)*(px)*(k) + (px)*(j+1) + (i);
    int pointer_jn = (py)*(px)*(k) + (px)*(j-1) + (i);
    int pointer_kp = (py)*(px)*(k+1) + (px)*(j) + (i);
    int pointer_kn = (py)*(px)*(k-1) + (px)*(j) + (i);

    double permitivity_ip=DDmaterial[pointer].k*DDmaterial[pointer_ip].k / (0.5*DDmaterial[pointer_ip].k+0.5*DDmaterial[pointer].k);
    double permitivity_in=DDmaterial[pointer].k*DDmaterial[pointer_in].k / (0.5*DDmaterial[pointer_in].k+0.5*DDmaterial[pointer].k);
    double permitivity_jp=DDmaterial[pointer].k*DDmaterial[pointer_jp].k / (0.5*DDmaterial[pointer_jp].k+0.5*DDmaterial[pointer].k);
    double permitivity_jn=DDmaterial[pointer].k*DDmaterial[pointer_jn].k / (0.5*DDmaterial[pointer_jn].k+0.5*DDmaterial[pointer].k);
    double permitivity_kp=DDmaterial[pointer].k*DDmaterial[pointer_kp].k / (0.5*DDmaterial[pointer_kp].k+0.5*DDmaterial[pointer].k);
    double permitivity_kn=DDmaterial[pointer].k*DDmaterial[pointer_kn].k / (0.5*DDmaterial[pointer_kn].k+0.5*DDmaterial[pointer].k);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
    double deltaz=abs(mesh[pointer_kp].coordZ-mesh[pointer_kn].coordZ)/2;
    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
    double zstep_p=abs(mesh[pointer_kp].coordZ-mesh[pointer].coordZ);
    double zstep_n=abs(mesh[pointer_kn].coordZ-mesh[pointer].coordZ);

    double f(0),df(0),phik(0);
    double volume=deltax*deltay*deltaz;

    phik=DDmaterial[pointer].phi;

    int Temp=DDmaterial[pointer].Type;

    switch(Temp){

    case 1: //1=channel
        f=volume*ni_nm*(DDmaterial[pointer].u*exp(phik/VT)-DDmaterial[pointer].v*exp((-1)*phik/VT)-DDmaterial[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(DDmaterial[pointer].u*exp(phik/VT)+DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 2: //2=insulator
        f=volume*DDmaterial[pointer].rho/e0;
        df=0;
        break;
    case 3: //3=electrolyte
        f=volume*C0*(DDmaterial[pointer].u*exp(phik/VT)-DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0;
       df=volume*C0*(DDmaterial[pointer].u*exp(phik/VT)+DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 4: //4=Substrate
        f=volume*ni_nm*(DDmaterial[pointer].u*exp(phik/VT)-DDmaterial[pointer].v*exp((-1)*phik/VT)-DDmaterial[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(DDmaterial[pointer].u*exp(phik/VT)+DDmaterial[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 5: //5=Analyte
        f=volume*DDmaterial[pointer].rho/e0;
        df=0;
        break;
    default:
        cout << "Error! Undefined Type for Poisson Solver @ PoissonGaussSeidelInner3D."<<endl;
        exit(0);
    }

    return ((permitivity_ip*DDmaterial[pointer_ip].phi/xstep_p+permitivity_in*DDmaterial[pointer_in].phi/xstep_n)*deltay*deltaz
           +(permitivity_jp*DDmaterial[pointer_jp].phi/ystep_p+permitivity_jn*DDmaterial[pointer_jn].phi/ystep_n)*deltax*deltaz
           +(permitivity_kp*DDmaterial[pointer_kp].phi/zstep_p+permitivity_kn*DDmaterial[pointer_kn].phi/zstep_n)*deltax*deltay + f - df*phik )
            /
            ((permitivity_ip/xstep_p+permitivity_in/xstep_n)*deltay*deltaz
            +(permitivity_jp/ystep_p+permitivity_jn/ystep_n)*deltax*deltaz
            +(permitivity_kp/zstep_p+permitivity_kn/zstep_n)*deltax*deltay - df);
}

double DDmodel::DD_ECSolver2D(){

    DD_loop=0;
    double errEC(0),errEC_max(0);

    do{
        DD_loop++;


        switch(StructureFlag){

        case 1:
            errEC=DD_ECTypeA2D();
            break;
        case 2:
            errEC=DD_ECTypeA2D();
            break;
        case 3:
            errEC=DD_ECTypeB2D();
            break;
        default:
            cout << "Not appropriate StructureFlag @ ECSolver2D." <<endl;
            exit(0);
        }

        if(errEC_max < errEC) {errEC_max=errEC;}

        if(DD_loop%1000==0)
        cout <<"EC:"<< DD_loop <<"\t" <<errEC <<"\t"<<errEC_max<<endl;

    }while(errEC>SimTolEC);

    return errEC_max;

}

double DDmodel::DD_ECTypeA2D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double uk = DDmaterial[pointer].u;

            DDmaterial[pointer].u=DD_ECInner2D(i,j);

            DDmaterial[pointer].r=DD_SRHrecomb2D( i,j);

            double error=VT*abs(log(DDmaterial[pointer].u)-log(uk));

            error=error/(abs(VT*log(uk))+1);

            if(error>max_val)
                max_val=error;
        }
    }

    DD_ECBC2D();

    return max_val;

}

double DDmodel::DD_ECTypeB2D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<SubstrateTOP; j++) {

            int pointer = (px)*(j) + (i);

            double uk = DDmaterial[pointer].u;

            DDmaterial[pointer].u=DD_ECInner2D(i,j);

            DDmaterial[pointer].r=DD_SRHrecomb2D( i,j);

            double error=abs(log(DDmaterial[pointer].u)-log(uk));

            if(error>max_val)
                max_val=error;
        }
    }

    /*
#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=ElectrolyteBottom+1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double uk = DDmaterial[pointer].u;

            DDmaterial[pointer].u=ECInner(i,j);

            DDmaterial[pointer].r=DD_SRHrecomb2D( i,j);

            double error=abs(log(DDmaterial[pointer].u)-log(uk));

            if(error>max_val)
                max_val=error;
        }
    }
*/
    DD_ECBC2D();

    return max_val;

}

double DDmodel::DD_ECInner2D(int i, int j){

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

    double volume;

    volume=deltax*deltay;

    double uf = DDmaterial[pointer].u;

    double Bip = (DDmaterial[pointer].mun+DDmaterial[pointer_ip].mun)/2.0*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT);
    double Bin = (DDmaterial[pointer].mun+DDmaterial[pointer_in].mun)/2.0*DD_Bern(DDmaterial[pointer_in].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_in].phi/VT);
    double Bjp = (DDmaterial[pointer].mun+DDmaterial[pointer_jp].mun)/2.0*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jp].phi/VT);
    double Bjn = (DDmaterial[pointer].mun+DDmaterial[pointer_jn].mun)/2.0*DD_Bern(DDmaterial[pointer_jn].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jn].phi/VT);

    double f=volume/VT*DDmaterial[pointer].r;

    double df=volume/VT*(pow(DDmaterial[pointer].v,2)*exp(-DDmaterial[pointer].phi/VT)+2*DDmaterial[pointer].v+exp(DDmaterial[pointer].phi/VT))/DDmaterial[pointer].tau/
                pow(DDmaterial[pointer].u*exp(DDmaterial[pointer].phi/VT)+DDmaterial[pointer].v*exp(-DDmaterial[pointer].phi/VT)+2,2);

    uf=((Bip*DDmaterial[pointer_ip].u/xstep_p+Bin*DDmaterial[pointer_in].u/xstep_n)*deltay
       +(Bjp*DDmaterial[pointer_jp].u/ystep_p+Bjn*DDmaterial[pointer_jn].u/ystep_n)*deltax
       -f+df*uf)/((Bip/xstep_p+Bin/xstep_n)*deltay+(Bjp/ystep_p+Bjn/ystep_n)*deltax+df);

    return uf;
}

double DDmodel::DD_HCSolver2D(){

    DD_loop=0;
    double errHC(0),errHC_max(0);

    do{
        DD_loop++;

        switch(StructureFlag){

        case 1:
            errHC=DD_HCTypeA2D();
            break;
        case 2:
            errHC=DD_HCTypeA2D();
            break;
        case 3:
            errHC=DD_HCTypeB2D();
            break;
        default:
            cout << "Not appropriate StructureFlag @ HCSolver2D." <<endl;
            exit(0);
        }

        if(errHC_max < errHC) {errHC_max=errHC;}

        if(DD_loop%1000==0)
        cout <<"HC:"<< DD_loop <<"\t" <<errHC <<"\t"<<errHC_max<<endl;

    }while(errHC>SimTolHC);

    return errHC_max;
}

double DDmodel::DD_HCTypeA2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double vk = DDmaterial[pointer].v;

            DDmaterial[pointer].v=DD_HCInner2D(i,j);

            DDmaterial[pointer].r=DD_SRHrecomb2D( i,j);

            double error=VT*abs(log(DDmaterial[pointer].v)-log(vk));

            error=error/(abs(VT*log(vk))+1);

            if(error>max_val)
                max_val=error;
        }
    }

    DD_HCBC2D();

    return max_val;
}

double DDmodel::DD_HCTypeB2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<SubstrateTOP; j++) {

            int pointer = (px)*(j) + (i);

            double vk = DDmaterial[pointer].v;

            DDmaterial[pointer].v=DD_HCInner2D(i,j);

            DDmaterial[pointer].r=DD_SRHrecomb2D( i,j);

            double error=abs(log(DDmaterial[pointer].v)-log(vk));

            if(error>max_val)
                max_val=error;
        }
    }

    /*
#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=ElectrolyteBottom+1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double vk = DDmaterial[pointer].v;

            DDmaterial[pointer].v=HCInner(i,j);

            DDmaterial[pointer].r=DD_SRHrecomb2D( i,j);

            double error=abs(log(DDmaterial[pointer].v)-log(vk));

            if(error>max_val)
                max_val=error;
        }
    }
    */
    DD_HCBC2D();

    return max_val;
}

double DDmodel::DD_HCInner2D(int i, int j){

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

    double volume;

    volume=deltax*deltay;

    double vf = DDmaterial[pointer].v;

    double Bip = (DDmaterial[pointer].mup+DDmaterial[pointer_ip].mup)/2.0*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);
    double Bin = (DDmaterial[pointer].mup+DDmaterial[pointer_in].mup)/2.0*DD_Bern(DDmaterial[pointer_in].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);
    double Bjp = (DDmaterial[pointer].mup+DDmaterial[pointer_jp].mup)/2.0*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);
    double Bjn = (DDmaterial[pointer].mup+DDmaterial[pointer_jn].mup)/2.0*DD_Bern(DDmaterial[pointer_jn].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);

    double f=volume/VT*DDmaterial[pointer].r;

    double df=volume/VT*(pow(DDmaterial[pointer].u,2)*exp(DDmaterial[pointer].phi/VT)+2*DDmaterial[pointer].u+exp(-DDmaterial[pointer].phi/VT))/DDmaterial[pointer].tau/
                pow(DDmaterial[pointer].u*exp(DDmaterial[pointer].phi/VT)+DDmaterial[pointer].v*exp(-DDmaterial[pointer].phi/VT)+2,2);

    vf =((Bip*DDmaterial[pointer_ip].v/xstep_p+Bin*DDmaterial[pointer_in].v/xstep_n)*deltay
         +(Bjp*DDmaterial[pointer_jp].v/ystep_p+Bjn*DDmaterial[pointer_jn].v/ystep_n)*deltax
         -f+df*vf)/((Bip/xstep_p+Bin/xstep_n)*deltay+(Bjp/ystep_p+Bjn/ystep_n)*deltax+df);

    return vf;
}

void DDmodel::DD_PoissonBC2D(){

    switch(StructureFlag){

    case 1:
        DD_PoissonBC2D_PN();
        break;
    case 2:
        DD_PoissonBC2D_MOSFET();
        break;
    case 3:
        DD_PoissonBC2D_ISFET();
        break;
    default:
        cout <<"Undefined Boundary @ PoissonBC2D."<<endl;
        exit(0);
    }
}

void DDmodel::DD_PoissonBC2D_PN(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);
            DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;

            pointer1 = (px)*(py-1) + (i);
            pointer2 = (px)*(py-2) + (i);
            DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
        }
}

void DDmodel::DD_PoissonBC2D_MOSFET(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);

            if( mesh[pointer1].coordX > JunctionLength && mesh[pointer1].coordX < lx-JunctionLength ){

                double ystep=abs(mesh[pointer1].coordY-mesh[pointer2].coordY);
                double qfactor=Si_permi/SiO2_permi*Tox/ystep;

                DDmaterial[pointer1].phi=(volG+qfactor*DDmaterial[pointer2].phi)/(1.0+qfactor);
            }
        }

#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);
            DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);
            DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
        }
}

void DDmodel::DD_PoissonBC2D_ISFET(){

    //substrate
    //if(mesh[pointer].coordY<=SubstrateThickness )

    //SD
    //if(mesh[pointer].coordY>SubstrateThickness-JunctionDepth && mesh[pointer].coordY<=SubstrateThickness )

    //oxide
    //if(mesh[pointer].coordY>SubstrateThickness && mesh[pointer].coordY<=SubstrateThickness+tox )

    //electrolyte
    //if(mesh[pointer].coordY>SubstrateThickness+tox)

#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;

            }

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                DDmaterial[pointer1].phi=DDmaterial[pointer2].phi;
            }
        }
}

void DDmodel::DD_ECBC2D(){

    switch(StructureFlag){

    case 1:
        DD_ECBC2D_PN();
        break;
    case 2:
        DD_ECBC2D_MOSFET();
        break;
    case 3:
        DD_ECBC2D_ISFET();
        break;
    default:
        cout <<"Undefined Boundary @ ECBC2D."<<endl;
        exit(0);
    }
}

void DDmodel::DD_HCBC2D(){

    switch(StructureFlag){

    case 1:
        DD_HCBC2D_PN();
        break;
    case 2:
        DD_HCBC2D_MOSFET();
        break;
    case 3:
        DD_HCBC2D_ISFET();
        break;
    default:
        cout <<"Undefined Boundary @ HCBC2D."<<endl;
        exit(0);
    }
}

void DDmodel::DD_ECBC2D_PN(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);
            DDmaterial[pointer1].u=DDmaterial[pointer2].u;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            pointer1 = (px)*(py-1) + (i);
            pointer2 = (px)*(py-2) + (i);
            DDmaterial[pointer1].u=DDmaterial[pointer2].u;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;
        }
}

void DDmodel::DD_HCBC2D_PN(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);
            DDmaterial[pointer1].v=DDmaterial[pointer2].v;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            pointer1 = (px)*(py-1) + (i);
            pointer2 = (px)*(py-2) + (i);
            DDmaterial[pointer1].v=DDmaterial[pointer2].v;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;
        }
}

void DDmodel::DD_ECBC2D_MOSFET(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);

            if( mesh[pointer1].coordX > JunctionLength && mesh[pointer1].coordX < lx-JunctionLength ){

                DDmaterial[pointer1].u=DDmaterial[pointer2].u;
                DDmaterial[pointer1].r=DDmaterial[pointer2].r;
            }
        }

#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);
            DDmaterial[pointer1].u=DDmaterial[pointer2].u;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);
            DDmaterial[pointer1].u=DDmaterial[pointer2].u;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;
        }
}

void DDmodel::DD_HCBC2D_MOSFET(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);

            if( mesh[pointer1].coordX > JunctionLength && mesh[pointer1].coordX < lx-JunctionLength ){

                DDmaterial[pointer1].v=DDmaterial[pointer2].v;
                DDmaterial[pointer1].r=DDmaterial[pointer2].r;
            }
        }

#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);
            DDmaterial[pointer1].v=DDmaterial[pointer2].v;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);
            DDmaterial[pointer1].v=DDmaterial[pointer2].v;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;
        }
}

void DDmodel::DD_ECBC2D_ISFET(){

    //substrate
    //if(mesh[pointer].coordY<=SubstrateThickness )

    //SD
    //if(mesh[pointer].coordY>SubstrateThickness-JunctionDepth && mesh[pointer].coordY<=SubstrateThickness )

    //oxide
    //if(mesh[pointer].coordY>SubstrateThickness && mesh[pointer].coordY<=SubstrateThickness+tox )

    //electrolyte
    //if(mesh[pointer].coordY>SubstrateThickness+tox)

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(SubstrateTOP) + (i);
            int pointer2 = (px)*(SubstrateTOP-1) + (i);
            DDmaterial[pointer1].u=DDmaterial[pointer2].u;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            pointer1 = (px)*(ElectrolyteBottom) + (i);
            pointer2 = (px)*(ElectrolyteBottom+1) + (i);
            DDmaterial[pointer1].u=DDmaterial[pointer2].u;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;
        }
#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                DDmaterial[pointer1].u=DDmaterial[pointer2].u;
                DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            }

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                DDmaterial[pointer1].u=DDmaterial[pointer2].u;
                DDmaterial[pointer1].r=DDmaterial[pointer2].r;
            }
        }
}

void DDmodel::DD_HCBC2D_ISFET(){

    //substrate
    //if(mesh[pointer].coordY<=SubstrateThickness )

    //SD
    //if(mesh[pointer].coordY>SubstrateThickness-JunctionDepth && mesh[pointer].coordY<=SubstrateThickness )

    //oxide
    //if(mesh[pointer].coordY>SubstrateThickness && mesh[pointer].coordY<=SubstrateThickness+tox )

    //electrolyte
    //if(mesh[pointer].coordY>SubstrateThickness+tox)

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(SubstrateTOP) + (i);
            int pointer2 = (px)*(SubstrateTOP-1) + (i);
            DDmaterial[pointer1].v=DDmaterial[pointer2].v;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            pointer1 = (px)*(ElectrolyteBottom) + (i);
            pointer2 = (px)*(ElectrolyteBottom+1) + (i);
            DDmaterial[pointer1].v=DDmaterial[pointer2].v;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;
        }
#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                DDmaterial[pointer1].v=DDmaterial[pointer2].v;
                DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            }

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                DDmaterial[pointer1].v=DDmaterial[pointer2].v;
                DDmaterial[pointer1].r=DDmaterial[pointer2].r;
            }
        }
}

void DDmodel::DD_ECBC3D(){

    switch(StructureFlag){

    case 1:
        DD_ECBC3D_PN();
        break;
    case 2:
        DD_ECBC3D_MOSFET();
        break;
    case 3:
        DD_ECBC3D_ISFET();
        break;
    case 4:
        DD_ECBC3D_1NWR();
        break;
    case 5:
        DD_ECBC3D_2NWR();
        break;
    default:
        cout <<"Undefined Boundary @ ECBC3D."<<endl;
        exit(0);
    }
}

void DDmodel::DD_ECBC3D_PN(){

    for(int i=0;i<px;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(0) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(1) +  (px)*(j) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(pz-1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(pz-2) +  (px)*(j) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(0) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(py-1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(py-2) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }
}

void DDmodel::DD_ECBC3D_MOSFET(){

    for(int i=0;i<px;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(0) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(1) +  (px)*(j) + (i);

            if(mesh[pointer0].coordX > JunctionLength && mesh[pointer0].coordX < lx-JunctionLength ){
                DDmaterial[pointer0].u=DDmaterial[pointer1].u;
                DDmaterial[pointer0].r=DDmaterial[pointer1].r;
            }
        }
    }

    for(int j=0;j<py;j++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(j) + (0);
            int pointer1 =(px)*(py)*(k) +  (px)*(j) + (1);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(pz-1) +  (px)*(j) + (px-1);
            pointer1 =(px)*(py)*(pz-2) +  (px)*(j) + (px-2);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(0) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(py-1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(py-2) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }
}

void DDmodel::DD_ECBC3D_ISFET(){

    for (int i=0; i<px; i++) {
        for (int j=0;j<py;j++){

            int pointer1 = (py)*(px)*(SubstrateTOP) + (px)*(j) + (i);
            int pointer2 = (py)*(px)*(SubstrateTOP-1) + (px)*(j) + (i);
            DDmaterial[pointer1].u=DDmaterial[pointer2].u;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            pointer1 = (py)*(px)*(ElectrolyteBottom) + (px)*(j) + (i);
            pointer2 = (py)*(px)*(ElectrolyteBottom+1) + (px)*(j) + (i);
            DDmaterial[pointer1].u=DDmaterial[pointer2].u;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;
        }
    }

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(mesh[pointer0].coordZ<=SubstrateThickness-JunctionDepth || mesh[pointer0].coordZ>SubstrateThickness ){
                DDmaterial[pointer0].u=DDmaterial[pointer1].u;
                DDmaterial[pointer0].r=DDmaterial[pointer1].r;
            }


            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(mesh[pointer0].coordZ<=SubstrateThickness-JunctionDepth || mesh[pointer0].coordZ>SubstrateThickness ){
                DDmaterial[pointer0].u=DDmaterial[pointer1].u;
                DDmaterial[pointer0].r=DDmaterial[pointer1].r;
            }
        }
    }
}

void DDmodel::DD_ECBC3D_1NWR(){

    for(int i=0;i<px;i++){
        for(int j=NWRleft1;j<=NWRright1;j++){
            int pointer0 =(px)*(py)*(NWRbottom1) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom1+1) +  (px)*(j) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop1-1) +  (px)*(j) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=NWRbottom1;k<=NWRtop1;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(NWRleft1) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(NWRleft1+1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(NWRright1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(NWRright1-1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }
}

void DDmodel::DD_ECBC3D_2NWR(){

    //NWR1
    for(int i=0;i<px;i++){
        for(int j=NWRleft1;j<=NWRright1;j++){
            int pointer0 =(px)*(py)*(NWRbottom1) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom1+1) +  (px)*(j) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop1-1) +  (px)*(j) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=NWRbottom1;k<=NWRtop1;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(NWRleft1) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(NWRleft1+1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(NWRright1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(NWRright1-1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    //NWR2
    for(int i=0;i<px;i++){
        for(int j=NWRleft2;j<=NWRright2;j++){
            int pointer0 =(px)*(py)*(NWRbottom2) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom2+1) +  (px)*(j) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop2) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop2-1) +  (px)*(j) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=NWRbottom2;k<=NWRtop2;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(NWRleft2) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(NWRleft2+1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(NWRright2) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(NWRright2-1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }
}

void DDmodel::DD_HCBC3D(){

    switch(StructureFlag){

    case 1:
        DD_HCBC3D_PN();
        break;
    case 2:
        DD_HCBC3D_MOSFET();
        break;
    case 3:
        DD_HCBC3D_ISFET();
        break;
    case 4:
        DD_HCBC3D_1NWR();
        break;
    case 5:
        DD_HCBC3D_2NWR();
        break;
    default:
        cout <<"Undefined Boundary @ ECBC2D."<<endl;
        exit(0);
    }
}

void DDmodel::DD_HCBC3D_PN(){

    for(int i=0;i<px;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(0) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(1) +  (px)*(j) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(pz-1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(pz-2) +  (px)*(j) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(0) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(1) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(py-1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(py-2) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }
}

void DDmodel::DD_HCBC3D_MOSFET(){

    for(int i=0;i<px;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(0) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(1) +  (px)*(j) + (i);

            if(mesh[pointer0].coordX > JunctionLength && mesh[pointer0].coordX < lx-JunctionLength ){
                DDmaterial[pointer0].v=DDmaterial[pointer1].v;
                DDmaterial[pointer0].r=DDmaterial[pointer1].r;
            }
        }
    }

    for(int j=0;j<py;j++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(j) + (0);
            int pointer1 =(px)*(py)*(k) +  (px)*(j) + (1);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(pz-1) +  (px)*(j) + (px-1);
            pointer1 =(px)*(py)*(pz-2) +  (px)*(j) + (px-2);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(0) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(1) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(py-1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(py-2) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }
}

void DDmodel::DD_HCBC3D_ISFET(){
    for (int i=0; i<px; i++) {
        for (int j=0;j<py;j++){

            int pointer1 = (py)*(px)*(SubstrateTOP) + (px)*(j) + (i);
            int pointer2 = (py)*(px)*(SubstrateTOP-1) + (px)*(j) + (i);
            DDmaterial[pointer1].v=DDmaterial[pointer2].v;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;

            pointer1 = (py)*(px)*(ElectrolyteBottom) + (px)*(j) + (i);
            pointer2 = (py)*(px)*(ElectrolyteBottom+1) + (px)*(j) + (i);
            DDmaterial[pointer1].v=DDmaterial[pointer2].v;
            DDmaterial[pointer1].r=DDmaterial[pointer2].r;
        }
    }

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            DDmaterial[pointer0].u=DDmaterial[pointer1].u;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(mesh[pointer1].coordZ<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordZ>SubstrateThickness){

                DDmaterial[pointer0].v=DDmaterial[pointer1].v;
                DDmaterial[pointer0].r=DDmaterial[pointer1].r;
            }


            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(mesh[pointer1].coordZ<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordZ>SubstrateThickness){

                DDmaterial[pointer0].v=DDmaterial[pointer1].v;
                DDmaterial[pointer0].r=DDmaterial[pointer1].r;
            }
        }
    }
}

void DDmodel::DD_HCBC3D_1NWR(){

    for(int i=0;i<px;i++){
        for(int j=NWRleft1;j<=NWRright1;j++){
            int pointer0 =(px)*(py)*(NWRbottom1) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom1+1) +  (px)*(j) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop1-1) +  (px)*(j) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=NWRbottom1;k<=NWRtop1;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(NWRleft1) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(NWRleft1+1) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(NWRright1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(NWRright1-1) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }
}

void DDmodel::DD_HCBC3D_2NWR(){

    //NWR1
    for(int i=NWRleft1;i<=NWRright1;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(NWRbottom1) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom1+1) +  (px)*(j) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop1-1) +  (px)*(j) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int k=NWRbottom1;k<=NWRtop1;k++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(k) +  (px)*(j) + (NWRleft1);
            int pointer1 =(px)*(py)*(k) +  (px)*(j) + (NWRleft1+1);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(j) + (NWRright1);
            pointer1 =(px)*(py)*(k) +  (px)*(j) + (NWRright1-1);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    //NWR2
    for(int i=NWRleft2;i<=NWRright2;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(NWRbottom2) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom2+1) +  (px)*(j) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop2) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop2-1) +  (px)*(j) + (i);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }

    for(int k=NWRbottom2;k<=NWRtop2;k++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(k) +  (px)*(j) + (NWRleft2);
            int pointer1 =(px)*(py)*(k) +  (px)*(j) + (NWRleft2+1);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(j) + (NWRright2);
            pointer1 =(px)*(py)*(k) +  (px)*(j) + (NWRright2-1);

            DDmaterial[pointer0].v=DDmaterial[pointer1].v;
            DDmaterial[pointer0].r=DDmaterial[pointer1].r;
        }
    }
}

void DDmodel::DD_PoissonBC3D(){

    switch(StructureFlag){

    case 1:
        DD_PoissonBC3D_PN();
        break;
    case 2:
        DD_PoissonBC3D_MOSFET();
        break;
    case 3:
        DD_PoissonBC3D_ISFET();
        break;
    case 4:
        DD_PoissonBC3D_1NWR();
        break;
    case 5:
        DD_PoissonBC3D_2NWR();
        break;
    default:
        cout <<"Undefined Boundary @ PoissonBC3D."<<endl;
        exit(0);
    }

}

void DDmodel::DD_PoissonBC3D_PN(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer0 = (py)*(px)*(0) + (px)*(j) + (i);
            int pointer1 = (py)*(px)*(1) + (px)*(j) + (i);
            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;

            pointer0 = (py)*(px)*(pz-1) + (px)*(j) + (i);
            pointer1 = (py)*(px)*(pz-2) + (px)*(j) + (i);
            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
        }
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);
            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);
            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
        }
    }
}

void DDmodel::DD_PoissonBC3D_MOSFET(){


#pragma omp parallel for
        for (int i=0; i<px; i++) {
            for (int j=0;j<py;j++){

                int pointer1 = (py)*(px)*(0) + (px)*(j) + (i);
                int pointer2 = (py)*(px)*(1) + (px)*(j) + (i);

                if(mesh[pointer1].coordX > JunctionLength && mesh[pointer1].coordX < lx-JunctionLength ){

                    double zstep=abs(mesh[pointer1].coordZ-mesh[pointer2].coordZ);
                    double qfactor=Si_permi/SiO2_permi*Tox/zstep;

                    DDmaterial[pointer1].phi=(volG+qfactor*DDmaterial[pointer2].phi)/(1.0+qfactor);
                }
            }
        }

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);
            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);
            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
        }
    }

    for (int k=0;k<pz;k++){
        for (int j=0;j<py;j++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);
            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);
            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
        }
    }
}

void DDmodel::DD_PoissonBC3D_ISFET(){

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(mesh[pointer0].coordZ<=SubstrateThickness-JunctionDepth || mesh[pointer0].coordZ>SubstrateThickness ){
                DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
            }

            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(mesh[pointer0].coordZ<=SubstrateThickness-JunctionDepth || mesh[pointer0].coordZ>SubstrateThickness ){
                DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
            }
        }
    }
}

void DDmodel::DD_PoissonBC3D_1NWR(){

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(DDmaterial[pointer0].Type!=1){
                DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
            }

            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(DDmaterial[pointer0].Type!=1){
                DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
            }
        }
    }
}

void DDmodel::DD_PoissonBC3D_2NWR(){

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(DDmaterial[pointer0].Type!=1){
                DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
            }

            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(DDmaterial[pointer0].Type!=1){
                DDmaterial[pointer0].phi=DDmaterial[pointer1].phi;
            }
        }
    }
}

void DDmodel::DD_BernoulliX(){

    double x(1.0),v1,v2;

    do{
        x=x/2.0;
        v1=x/(exp(x)-1.0);
        v2=1.0-x/2.0;

    }while((v1-v2)>1e-10);
    bernXl=x;
}

void DDmodel::DD_PrintMaterial2D(string path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tY(2)\tK(3)\tdop(4)\tphi(5)\tu(6)\tv(7)\tr(8)\trho(9)\tmun(10)\tmup(11)\ttau(12)\tEx(13)\tEy(14)\tType(15)\tCrho(16)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            output << mesh[pointer].coordX << '\t' << mesh[pointer].coordY << '\t'
                   << DDmaterial[pointer].k << '\t' <<DDmaterial[pointer].dop << '\t' <<DDmaterial[pointer].phi << '\t'
                   << DDmaterial[pointer].u << '\t' << DDmaterial[pointer].v << '\t' << DDmaterial[pointer].r << '\t'
                   << DDmaterial[pointer].rho << '\t'<< DDmaterial[pointer].mun << '\t'<< DDmaterial[pointer].mup << '\t'
                   << DDmaterial[pointer].tau << '\t'<< DDmaterial[pointer].Ex << '\t'<< DDmaterial[pointer].Ey << '\t'
                   << DDmaterial[pointer].Type << '\t'<< DDmaterial[pointer].Crho <<endl;

        }
    }

    output.close();
}

void DDmodel::DD_PrintMaterial3D(string path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tY(2)\tZ(3)\tK(4)\tdop(5)\tphi(6)\tu(7)\tv(8)\tr(9)\trho(10)\tmun(11)\tmup(12)\ttau(13)\tEx(14)\tEy(15)\tEz(16)\tType(17)\tCrho(18)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){
                int pointer = (px)*(py)*(k) + (px)*(j) + (i);
                output << mesh[pointer].coordX << '\t' << mesh[pointer].coordY << '\t'<< mesh[pointer].coordZ << '\t'
                       << DDmaterial[pointer].k << '\t' <<DDmaterial[pointer].dop << '\t' <<DDmaterial[pointer].phi << '\t'
                       << DDmaterial[pointer].u << '\t' << DDmaterial[pointer].v << '\t' << DDmaterial[pointer].r << '\t'
                       << DDmaterial[pointer].rho << '\t'<< DDmaterial[pointer].mun << '\t'<< DDmaterial[pointer].mup << '\t'
                       << DDmaterial[pointer].tau << '\t'<< DDmaterial[pointer].Ex << '\t'<< DDmaterial[pointer].Ey << '\t'
                       << DDmaterial[pointer].Ez << '\t'<< DDmaterial[pointer].Type << '\t'<< DDmaterial[pointer].Crho <<endl;
            }
        }
    }

    output.close();
}

void DDmodel::DD_ReadMaterial2D(string path){

    DD_Initialize();

    fstream input;
    input.open (path, fstream::in);
    if(!input){
        cout << "File Reading Fail. File name : "<<path<<endl;
    }

    for (int i=0;i<3;i++) input.ignore(256,'#');

    double buffer;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            input >> buffer >> buffer
                  >> DDmaterial[pointer].k >> DDmaterial[pointer].dop >> DDmaterial[pointer].phi
                  >> DDmaterial[pointer].u >> DDmaterial[pointer].v >> DDmaterial[pointer].r
                  >> DDmaterial[pointer].rho >> DDmaterial[pointer].mun >> DDmaterial[pointer].mup
                  >> DDmaterial[pointer].tau >> DDmaterial[pointer].Ex >> DDmaterial[pointer].Ey
                  >> DDmaterial[pointer].Type >> DDmaterial[pointer].Crho;
        }
    }

    input.close();
}

void DDmodel::DD_ReadMaterial3D(string path){

    DD_Initialize();

    fstream input;
    input.open (path, fstream::in);
    if(!input){
        cout << "File Reading Fail. File name : "<<path<<endl;
    }

    for (int i=0;i<3;i++) input.ignore(256,'#');

    double buffer;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            input >> buffer >> buffer >> buffer
                  >> DDmaterial[pointer].k >> DDmaterial[pointer].dop >> DDmaterial[pointer].phi
                  >> DDmaterial[pointer].u >> DDmaterial[pointer].v >> DDmaterial[pointer].r
                  >> DDmaterial[pointer].rho >> DDmaterial[pointer].mun >> DDmaterial[pointer].mup
                  >> DDmaterial[pointer].tau >> DDmaterial[pointer].Ex >> DDmaterial[pointer].Ey
                  >> DDmaterial[pointer].Ez >> DDmaterial[pointer].Type >> DDmaterial[pointer].Crho;
        }
    }

    input.close();
}

void DDmodel::DD_RhoCalculation2D(){

    #pragma omp parallel for
    for (int i=1;i<px-1;i++){
        for (int j=1;j<py-1;j++){
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

            double permitivity_ip=DDmaterial[pointer].k*DDmaterial[pointer_ip].k / (0.5*DDmaterial[pointer_ip].k+0.5*DDmaterial[pointer].k);
            double permitivity_in=DDmaterial[pointer].k*DDmaterial[pointer_in].k / (0.5*DDmaterial[pointer_in].k+0.5*DDmaterial[pointer].k);
            double permitivity_jp=DDmaterial[pointer].k*DDmaterial[pointer_jp].k / (0.5*DDmaterial[pointer_jp].k+0.5*DDmaterial[pointer].k);
            double permitivity_jn=DDmaterial[pointer].k*DDmaterial[pointer_jn].k / (0.5*DDmaterial[pointer_jn].k+0.5*DDmaterial[pointer].k);

            double charge_totalemp = (-1)*((permitivity_ip*DDmaterial[pointer_ip].phi - permitivity_ip*DDmaterial[pointer].phi)/xstep_p*deltay
                                          +(permitivity_in*DDmaterial[pointer_in].phi - permitivity_in*DDmaterial[pointer].phi)/xstep_n*deltay
                                          +(permitivity_jp*DDmaterial[pointer_jp].phi - permitivity_jp*DDmaterial[pointer].phi)/ystep_p*deltax
                                          +(permitivity_jn*DDmaterial[pointer_jn].phi - permitivity_jn*DDmaterial[pointer].phi)/ystep_n*deltax)
                                          /(deltax*deltay);

            DDmaterial[pointer].Crho=charge_totalemp*e0;
        }
    }
}

void DDmodel::DD_RhoCalculation3D(){

    #pragma omp parallel for
    for (int i=1;i<px-1;i++){
        for (int j=1;j<py-1;j++){
            for (int k=1;k<pz-1;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);
                int pointer_ip = (py)*(px)*(k) + (px)*(j) + (i+1);
                int pointer_in = (py)*(px)*(k) + (px)*(j) + (i-1);
                int pointer_jp = (py)*(px)*(k) + (px)*(j+1) + (i);
                int pointer_jn = (py)*(px)*(k) + (px)*(j-1) + (i);
                int pointer_kp = (py)*(px)*(k+1) + (px)*(j) + (i);
                int pointer_kn = (py)*(px)*(k-1) + (px)*(j) + (i);

                double permitivity_ip=DDmaterial[pointer].k*DDmaterial[pointer_ip].k / (0.5*DDmaterial[pointer_ip].k+0.5*DDmaterial[pointer].k);
                double permitivity_in=DDmaterial[pointer].k*DDmaterial[pointer_in].k / (0.5*DDmaterial[pointer_in].k+0.5*DDmaterial[pointer].k);
                double permitivity_jp=DDmaterial[pointer].k*DDmaterial[pointer_jp].k / (0.5*DDmaterial[pointer_jp].k+0.5*DDmaterial[pointer].k);
                double permitivity_jn=DDmaterial[pointer].k*DDmaterial[pointer_jn].k / (0.5*DDmaterial[pointer_jn].k+0.5*DDmaterial[pointer].k);
                double permitivity_kp=DDmaterial[pointer].k*DDmaterial[pointer_kp].k / (0.5*DDmaterial[pointer_kp].k+0.5*DDmaterial[pointer].k);
                double permitivity_kn=DDmaterial[pointer].k*DDmaterial[pointer_kn].k / (0.5*DDmaterial[pointer_kn].k+0.5*DDmaterial[pointer].k);

                double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
                double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
                double deltaz=abs(mesh[pointer_kp].coordZ-mesh[pointer_kn].coordZ)/2;
                double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
                double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
                double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                double zstep_p=abs(mesh[pointer_kp].coordZ-mesh[pointer].coordZ);
                double zstep_n=abs(mesh[pointer_kn].coordZ-mesh[pointer].coordZ);

                double charge_totalemp = (-1)*((permitivity_ip*DDmaterial[pointer_ip].phi - permitivity_ip*DDmaterial[pointer].phi)/xstep_p*deltay*deltaz
                                              +(permitivity_in*DDmaterial[pointer_in].phi - permitivity_in*DDmaterial[pointer].phi)/xstep_n*deltay*deltaz
                                              +(permitivity_jp*DDmaterial[pointer_jp].phi - permitivity_jp*DDmaterial[pointer].phi)/ystep_p*deltax*deltaz
                                              +(permitivity_jn*DDmaterial[pointer_jn].phi - permitivity_jn*DDmaterial[pointer].phi)/ystep_n*deltax*deltaz
                                              +(permitivity_kp*DDmaterial[pointer_kp].phi - permitivity_kp*DDmaterial[pointer].phi)/zstep_p*deltax*deltay
                                              +(permitivity_kn*DDmaterial[pointer_kn].phi - permitivity_kn*DDmaterial[pointer].phi)/zstep_n*deltax*deltay)
                                              /(deltax*deltay*deltaz);

                DDmaterial[pointer].Crho=charge_totalemp*e0;
                //Crho = [C/nm3]
            }
        }
    }
}

void DDmodel::DD_EfieldCalculation2D(){

    #pragma omp parallel for
    for (int i=1;i<px-1;i++){
        for (int j=1;j<py-1;j++){

            int pointer = (px)*(j) + (i);
            int pointer_ip =   (px)*(j) + (i+1);
            int pointer_in =   (px)*(j) + (i-1);
            int pointer_jp =   (px)*(j+1) + (i);
            int pointer_jn =   (px)*(j-1) + (i);

            double Efx=(-1)*(DDmaterial[pointer_ip].phi-DDmaterial[pointer_in].phi)/((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)*1e-9);
            double Efy=(-1)*(DDmaterial[pointer_jp].phi-DDmaterial[pointer_jn].phi)/((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)*1e-9);

            DDmaterial[pointer].Ex=Efx;
            DDmaterial[pointer].Ey=Efy;
        }
    }

    for (int i=0;i<px;i++){

        int pointer1 = (px)*(0) + (i);
        int pointer2 = (px)*(1) + (i);
        DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
        DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;

        pointer1 = (px)*(py-1) + (i);
        pointer2 = (px)*(py-2) + (i);
        DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
        DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;
    }

    for (int j=0;j<py;j++){

        int pointer1 = (px)*(j) + (0);
        int pointer2 = (px)*(j) + (1);
        DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
        DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;

        pointer1 = (px)*(j) + (px-1);
        pointer2 = (px)*(j) + (px-2);
        DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
        DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;
    }

}

void DDmodel::DD_EfieldCalculation3D(){

    #pragma omp parallel for
    for (int i=1;i<px-1;i++){
        for (int j=1;j<py-1;j++){
            for (int k=1;k<pz-1;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);
                int pointer_ip = (py)*(px)*(k) + (px)*(j) + (i+1);
                int pointer_in = (py)*(px)*(k) + (px)*(j) + (i-1);
                int pointer_jp = (py)*(px)*(k) + (px)*(j+1) + (i);
                int pointer_jn = (py)*(px)*(k) + (px)*(j-1) + (i);
                int pointer_kp = (py)*(px)*(k+1) + (px)*(j) + (i);
                int pointer_kn = (py)*(px)*(k-1) + (px)*(j) + (i);

                double Efx=(-1)*(DDmaterial[pointer_ip].phi-DDmaterial[pointer_in].phi)/((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)*1e-9);
                double Efy=(-1)*(DDmaterial[pointer_jp].phi-DDmaterial[pointer_jn].phi)/((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)*1e-9);
                double Efz=(-1)*(DDmaterial[pointer_kp].phi-DDmaterial[pointer_kn].phi)/((mesh[pointer_kp].coordZ-mesh[pointer_kn].coordZ)*1e-9);

                DDmaterial[pointer].Ex=Efx;
                DDmaterial[pointer].Ey=Efy;
                DDmaterial[pointer].Ez=Efz;
            }
        }
    }

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer1 = (py)*(px)*(0) + (px)*(j) + (i);
            int pointer2 = (py)*(px)*(1) + (px)*(j) + (i);
            DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
            DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;
            DDmaterial[pointer1].Ez=DDmaterial[pointer2].Ez;

            pointer1 = (py)*(px)*(pz-1) + (px)*(j) + (i);
            pointer2 = (py)*(px)*(pz-2) + (px)*(j) + (i);
            DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
            DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;
            DDmaterial[pointer1].Ez=DDmaterial[pointer2].Ez;
        }
    }

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer1 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer2 = (py)*(px)*(k) + (px)*(1) + (i);
            DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
            DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;
            DDmaterial[pointer1].Ez=DDmaterial[pointer2].Ez;

            pointer1 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer2 = (py)*(px)*(k) + (px)*(py-2) + (i);
            DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
            DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;
            DDmaterial[pointer1].Ez=DDmaterial[pointer2].Ez;
        }
    }

    for (int k=0;k<pz;k++){
        for (int j=0;j<py;j++){

            int pointer1 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer2 = (py)*(px)*(k) + (px)*(j) + (1);
            DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
            DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;
            DDmaterial[pointer1].Ez=DDmaterial[pointer2].Ez;

            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer2 = (py)*(px)*(k) + (px)*(j) + (px-2);
            DDmaterial[pointer1].Ex=DDmaterial[pointer2].Ex;
            DDmaterial[pointer1].Ey=DDmaterial[pointer2].Ey;
            DDmaterial[pointer1].Ez=DDmaterial[pointer2].Ez;
        }
    }
}

void DDmodel::DD_Initialize(){

    #pragma omp parallel for
    for(int i=0;i<L;i++){
        DDmaterial[i].Crho=0;
        DDmaterial[i].dop=0;
        DDmaterial[i].Ex=0;
        DDmaterial[i].Ey=0;
        DDmaterial[i].Ez=0;
        DDmaterial[i].k=0;
        DDmaterial[i].mun=0;
        DDmaterial[i].mup=0;
        DDmaterial[i].phi=0;
        DDmaterial[i].r=0;
        DDmaterial[i].rho=0;
        DDmaterial[i].tau=0;
        DDmaterial[i].Type=0;
        DDmaterial[i].u=0;
        DDmaterial[i].v=0;
    }
}

double DDmodel::DD_munCal(double T, double dopping, int f){

    // 1200 cm^2/v.s  0.12 m^2/v.s  0.12*1e18 nm^2/v.s

    if(f==1){ //channel
        return (68.5+(1414-68.5)/(1+pow(dopping/9.2e16,0.711)))*1e14;
    }
    else if(f==2){ //oxide
        return 0;
    }
    else if(f==3){ //electrolyte
        return (68.5+(1414-68.5)/(1+pow(dopping/9.2e16,0.711)))*1e14;
    }
    else if(f==4){ //sub
        return (68.5+(1414-68.5)/(1+pow(dopping/9.2e16,0.711)))*1e14;
    }
    else{
        cout << "undefined behavior @ munCal."<<endl;
        exit(0);
    }
    //return  68.5 + (1414-68.5)/(1+pow(dopping/9.2e16,0.711))
    //return (88*pow(T/300,-0.57)+1252*pow(T/300,-2.33)/(1+N/(1.432e23*pow(T/300,2.546))))*1e-4;
}

double DDmodel::DD_mupCal(double T,double dopping, int f){

    // 1200 cm^2/v.s  0.12 m^2/v.s  0.12*1e18 nm^2/v.s

    if(f==1){ //channel
        return (44.9+(470.5-44.9)/(1+pow(dopping/2.23e17,0.719)))*1e14;
    }
    else if(f==2){ //oxide
        return 0;
    }
    else if(f==3){ //electrolyte
        return (44.9+(470.5-44.9)/(1+pow(dopping/2.23e17,0.719)))*1e14;
    }
    else if(f==4){ //sub
        return (44.9+(470.5-44.9)/(1+pow(dopping/2.23e17,0.719)))*1e14;
    }
    else{
        cout << "undefined behavior @ mupCal."<<endl;
        exit(0);
    }
    //return (54.3*pow(T/300,-0.57)+407*pow(T/300,-2.33)/(1+N/(2.67e23*pow(T/300,2.546))))*1e-4
}

double DDmodel::DD_tauNCal(double dopping){

    //Nd dopping
    //calculate hole lifetime
    return 3.94e-4/(1+dopping/1e21);
}

double DDmodel::DD_tauPCal(double dopping){

    //Na dopping
    //calculate electron lifetime
    return 3.94e-4/(1+dopping/1e21);
}

double DDmodel::DD_SRHrecomb2D(int i, int j){

    int pointer = (px)*(j) + (i);
    //http://www.iue.tuwien.ac.at/phd/entner/node11.html
    return (DDmaterial[pointer].u*DDmaterial[pointer].v-1)/DDmaterial[pointer].tau/(DDmaterial[pointer].u*exp(DDmaterial[pointer].phi/VT)+DDmaterial[pointer].v*exp((-1)*DDmaterial[pointer].phi/VT)+2);
}

double DDmodel::DD_SRHrecomb3D(int i, int j, int k){

    int pointer = (px)*(py)*(k) + (px)*(j) + (i);
    //http://www.iue.tuwien.ac.at/phd/entner/node11.html
    return (DDmaterial[pointer].u*DDmaterial[pointer].v-1)/DDmaterial[pointer].tau/(DDmaterial[pointer].u*exp(DDmaterial[pointer].phi/VT)+DDmaterial[pointer].v*exp((-1)*DDmaterial[pointer].phi/VT)+2);
}

double DDmodel::DD_Bern(double dphi, double phi)
{
    if(abs(dphi)<bernXl){
        return exp(phi)*(1.0-dphi/2);
    }
    else{
        return exp(phi)*dphi/(exp(dphi)-1.0);
    }
}

void DDmodel::DD_IdVG2D(){

    int numIter(0);
    double errMax(0),errPhi(0),errElec(0),errHole(0);
    int iter_Phi(0),iter_Elec(0),iter_Hole(0);
    int index(0);
    ofstream  output1;
    ofstream  output2;

    output1.open("current.txt", fstream::out | fstream::trunc);
    output1.precision(6);
    output1<<"VoltS="<<volS<<"\t"<<"VoltD="<<volD<<"\t"<<"SimTolPoisson="<<SimTolPoisson<<"\t"<<endl;
    output1<<"Vs(1)"<<"\t"<<"Vg(2)"<<"\t"<<"Vd(3)"<<"\t"<<"J_Sn(A/nm)(4)"<<"\t"<<"J_Sp(A/nm)(5)"<<"\t"<<"J_Dn(A/nm)(6)"<<"\t"<<"J_Dp(A/nm)(7)"
           <<"\t"<<"J_S(A/nm)(8)"<<"\t"<<"J_D(A/nm)(9)"<<"\t"<<"J_Bn(A/nm)(10)"<<"\t"<<"J_Bp(A/nm)(11)"<<"\t"<<"J_B(A/nm)(12)"<<endl;
    output1<<"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="<<endl;
    output1.close();


    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);

    do{
        volG=(volGi+index*volGs);
        DD_InitialGuess2D();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;

        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=DD_PoissonSolver2D();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=DD_ECSolver2D();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<"\t";

            //hole===============
            errHole=DD_HCSolver2D();
            iter_Hole=DD_loop;
            if(errHole>errMax)
                errMax=errHole;

            output2<<"Hole:" << iter_Hole <<"\t"<<errHole<<"\t";
            output2<<"Max.err:" <<errMax <<endl;

        }while( (iter_Hole!=1 || iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        DD_RhoCalculation2D();
        DD_EfieldCalculation2D();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<".txt";
        name2=name1.str();
        DD_PrintMaterial2D(name2.c_str());

        DD_Jcal2D();

        index++;
        numIter=0;

    }while(volGi+index*volGs<(volGe+0.001));

    output2.close();
    cout << "Simulation Process Finished."<<endl;

}

void DDmodel::DD_IdVD2D(){

    int numIter(0);
    double errMax(0),errPhi(0),errElec(0),errHole(0);
    int iter_Phi(0),iter_Elec(0),iter_Hole(0);
    int index(0);
    ofstream  output1;
    ofstream  output2;

    output1.open("current.txt", fstream::out | fstream::trunc);
    output1.precision(6);
    output1<<"Vs(1)"<<"\t"<<"Vg(2)"<<"\t"<<"Vd(3)"<<"\t"<<"J_Sn(A/nm)(4)"<<"\t"<<"J_Sp(A/nm)(5)"<<"\t"<<"J_Dn(A/nm)(6)"<<"\t"<<"J_Dp(A/nm)(7)"
           <<"\t"<<"J_S(A/nm)(8)"<<"\t"<<"J_D(A/nm)(9)"<<"\t"<<"J_Bn(A/nm)(10)"<<"\t"<<"J_Bp(A/nm)(11)"<<"\t"<<"J_B(A/nm)(12)"<<endl;
    output1<<"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="<<endl;
    output1.close();

    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);

    do{
        volD=(volDi+index*volDs);
        DD_InitialGuess2D();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;


        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=DD_PoissonSolver2D();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=DD_ECSolver2D();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<"\t";

            //hole===============
            errHole=DD_HCSolver2D();
            iter_Hole=DD_loop;
            if(errHole>errMax)
                errMax=errHole;

            output2<<"Hole:" << iter_Hole <<"\t"<<errHole<<"\t";
            output2<<"Max.err:" <<errMax <<endl;

        }while( (iter_Hole!=1 || iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        DD_RhoCalculation2D();
        DD_EfieldCalculation2D();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<".txt";
        name2=name1.str();
        DD_PrintMaterial2D(name2.c_str());

        DD_Jcal2D();

        index++;
        numIter=0;

    }while(volDi+index*volDs<volDe+0.001);

    output2.close();
    cout << "Simulation Process Finished."<<endl;
}

void DDmodel::DD_IdVG3D(){

    int numIter(0);
    double errMax(0),errPhi(0),errElec(0),errHole(0);
    int iter_Phi(0),iter_Elec(0),iter_Hole(0);
    int index(0);
    ofstream  output1;
    ofstream  output2;

    output1.open("current.txt", fstream::out | fstream::trunc);
    output1.precision(6);
    output1<<"VoltS="<<volS<<"\t"<<"VoltD="<<volD<<"\t"<<"SimTolPoisson="<<SimTolPoisson<<"\t"<<endl;
    output1<<"Vs(1)"<<"\t"<<"Vg(2)"<<"\t"<<"Vd(3)"<<"\t"<<"J_Sn(A/nm)(4)"<<"\t"<<"J_Sp(A/nm)(5)"<<"\t"<<"J_Dn(A/nm)(6)"<<"\t"<<"J_Dp(A/nm)(7)"
           <<"\t"<<"J_S(A/nm)(8)"<<"\t"<<"J_D(A/nm)(9)"<<"\t"<<"J_Bn(A/nm)(10)"<<"\t"<<"J_Bp(A/nm)(11)"<<"\t"<<"J_B(A/nm)(12)"<<endl;
    output1<<"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="<<endl;
    output1.close();


    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);

    do{
        volG=(volGi+index*volGs);
        DD_InitialGuess3D();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;

        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=DD_PoissonSolver3D();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=DD_ECSolver3D();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<"\t";

            //hole===============
            errHole=DD_HCSolver3D();
            iter_Hole=DD_loop;
            if(errHole>errMax)
                errMax=errHole;

            output2<<"Hole:" << iter_Hole <<"\t"<<errHole<<"\t";
            output2<<"Max.err:" <<errMax <<endl;

        }while( (iter_Hole!=1 || iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        DD_RhoCalculation3D();
        DD_EfieldCalculation3D();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<".txt";
        name2=name1.str();
        DD_PrintMaterial3D(name2.c_str());

        DD_Jcal3D();

        index++;
        numIter=0;

    }while(volGi+index*volGs<(volGe+0.001));

    output2.close();
    cout << "Simulation Process Finished."<<endl;

}

void DDmodel::DD_IdVD3D(){

    int numIter(0);
    double errMax(0),errPhi(0),errElec(0),errHole(0);
    int iter_Phi(0),iter_Elec(0),iter_Hole(0);
    int index(0);
    ofstream  output1;
    ofstream  output2;

    output1.open("current.txt", fstream::out | fstream::trunc);
    output1.precision(6);
    output1<<"Vs(1)"<<"\t"<<"Vg(2)"<<"\t"<<"Vd(3)"<<"\t"<<"J_Sn(A/nm)(4)"<<"\t"<<"J_Sp(A/nm)(5)"<<"\t"<<"J_Dn(A/nm)(6)"<<"\t"<<"J_Dp(A/nm)(7)"
           <<"\t"<<"J_S(A/nm)(8)"<<"\t"<<"J_D(A/nm)(9)"<<"\t"<<"J_Bn(A/nm)(10)"<<"\t"<<"J_Bp(A/nm)(11)"<<"\t"<<"J_B(A/nm)(12)"<<endl;
    output1<<"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="<<endl;
    output1.close();

    output2.open("convergence.txt", fstream::out | fstream::trunc);
    output2.precision(10);

    do{
        volD=(volDi+index*volDs);
        DD_InitialGuess3D();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;


        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=DD_PoissonSolver3D();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=DD_ECSolver3D();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<"\t";

            //hole===============
            errHole=DD_HCSolver3D();
            iter_Hole=DD_loop;
            if(errHole>errMax)
                errMax=errHole;

            output2<<"Hole:" << iter_Hole <<"\t"<<errHole<<"\t";
            output2<<"Max.err:" <<errMax <<endl;

        }while( (iter_Hole!=1 || iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        DD_RhoCalculation3D();
        DD_EfieldCalculation3D();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<".txt";
        name2=name1.str();
        DD_PrintMaterial3D(name2.c_str());

        DD_Jcal3D();

        index++;
        numIter=0;

    }while(volDi+index*volDs<volDe+0.001);

    output2.close();
    cout << "Simulation Process Finished."<<endl;
}


void DDmodel::DD_Jcal2D(){

    switch(StructureFlag){

    case 1:
        DD_Jcal2D_PN();
        break;
    case 2:
        DD_Jcal2D_MOSFET();
        break;
    case 3:
        DD_Jcal2D_ISFET();
        break;
    default:
        cout << "Not appropriate StructureFlag @ DD_Jcal2D."<<endl;
        exit(0);
    }
}

void DDmodel::DD_Jcal2D_PN(void){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0);

    DD_JcalS2D_PN(Current_Sn,Current_Sp);
    DD_JcalD2D_PN(Current_Dn,Current_Dp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific<<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp<<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<endl;
    output1.close();
}

void DDmodel::DD_JcalS2D_PN(double &JSn,double &JSp){

    for (int j=0; j<py; j++) {
        int pointer = (px)*(j) + (0);
        int pointer_ip = (px)*(j) + (1);
        int pointer_jp = (px)*(j+1) + (0);
        int pointer_jn = (px)*(j-1) + (1);

        double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);

        double deltay=0;

        if(j==0){
            deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
        }
        else if(j==py-1){
            deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
        }else {
            deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
        }

        JSn+= (DDmaterial[pointer].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
             (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay/xstep;

        JSp+= (DDmaterial[pointer].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
             (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay/xstep;
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalD2D_PN(double &JDn, double &JDp){

    for (int j=0; j<py; j++) {
        int pointer = (px)*(j) + (px-2);
        int pointer_ip = (px)*(j) + (px-1);
        int pointer_jp = (px)*(j+1) + (px-2);
        int pointer_jn = (px)*(j-1) + (px-1);

        double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);

        double deltay=0;

        if(j==0){
            deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
        }
        else if(j==py-1){
            deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
        }else {
            deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
        }

        JDn+= (DDmaterial[pointer].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
             (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay/xstep;

        JDp+= (DDmaterial[pointer].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
             (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay/xstep;
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::DD_Jcal2D_MOSFET(void){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    DD_JcalS2D_MOSFET(Current_Sn,Current_Sp);
    DD_JcalD2D_MOSFET(Current_Dn,Current_Dp);
    DD_JcalB2D_MOSFET(Current_Bn,Current_Bp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific
          <<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
          <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"
          <<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;
    output1.close();
}

void DDmodel::DD_JcalS2D_MOSFET(double &JSn, double &JSp){

    for (int i=0; i<px; i++) {
        int pointer = (px)*(0) + (i);

        if(mesh[pointer].coordX<=JunctionLength){
            int pointer_jp = (px)*(1) + (i);
            int pointer_ip = (px)*(0) + (i+1);
            int pointer_in = (px)*(0) + (i-1);
            double ystep=abs(mesh[pointer].coordY-mesh[pointer_jp].coordY);
            double deltax=0;

            if(i==0){
                deltax=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
            }
            else if(i==px-1){
                deltax=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
            }else{
                deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
            }

            JSn+= (DDmaterial[pointer_jp].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jp].phi/VT)*
                 (DDmaterial[pointer_jp].u-DDmaterial[pointer].u)*deltax/ystep;

            JSp+= (DDmaterial[pointer_jp].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                 (DDmaterial[pointer_jp].v-DDmaterial[pointer].v)*deltax/ystep;
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalD2D_MOSFET(double &JDn, double &JDp){

    for (int i=0; i<px; i++) {
        int pointer = (px)*(0) + (i);

        if(mesh[pointer].coordX>=lx-JunctionLength){
            int pointer_jp = (px)*(1) + (i);
            int pointer_ip = (px)*(0) + (i+1);
            int pointer_in = (px)*(0) + (i-1);
            double ystep=abs(mesh[pointer].coordY-mesh[pointer_jp].coordY);
            double deltax=0;

            if(i==0){
                deltax=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
            }
            else if(i==px-1){
                deltax=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
            }else{
                deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
            }

            JDn+= (DDmaterial[pointer_jp].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jp].phi/VT)*
                 (DDmaterial[pointer_jp].u-DDmaterial[pointer].u)*deltax/ystep;

            JDp+= (DDmaterial[pointer_jp].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                 (DDmaterial[pointer_jp].v-DDmaterial[pointer].v)*deltax/ystep;
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalB2D_MOSFET(double &JBn, double &JBp){

    for (int i=0; i<px; i++) {
        int pointer = (px)*(py-2) + (i);
        int pointer_jp = (px)*(py-1) + (i);
        int pointer_ip = (px)*(py-2) + (i+1);
        int pointer_in = (px)*(py-2) + (i-1);

        double ystep=abs(mesh[pointer].coordY-mesh[pointer_jp].coordY);

        double deltax=0;

        if(i==0){
            deltax=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
        }
        else if(i==px-1){
            deltax=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
        }else{
            deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
        }

        JBn+= (DDmaterial[pointer_jp].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jp].phi/VT)*
             (DDmaterial[pointer_jp].u-DDmaterial[pointer].u)*deltax/ystep;

        JBp+= (DDmaterial[pointer].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
             (DDmaterial[pointer_jp].v-DDmaterial[pointer].v)*deltax/ystep;
    }

    JBn=JBn*q0*ni_nm*VT*(-1);
    JBp=JBp*q0*ni_nm*VT;
}

void DDmodel::DD_Jcal2D_ISFET(void){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    DD_JcalS2D_ISFET(Current_Sn,Current_Sp);
    DD_JcalD2D_ISFET(Current_Dn,Current_Dp);
    DD_JcalB2D_ISFET(Current_Bn,Current_Bp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific
          <<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
          <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"
          <<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;
    output1.close();
}

void DDmodel::DD_JcalS2D_ISFET(double &JSn, double &JSp){

    for (int j=0; j<py; j++) {
        int pointer = (px)*(j) + (0);
        int pointer_ip = (px)*(j) + (1);
        int pointer_jp = (px)*(j+1) + (0);
        int pointer_jn = (px)*(j-1) + (1);

        double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);

        double deltay=0;
        double ystep_p=0;
        double ystep_n=0;


        if(mesh[pointer].coordY>SubstrateThickness-JunctionLength && mesh[pointer].coordY<=SubstrateThickness){

            if(j==0){
                ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                deltay=ystep_p;
            }
            else if(j==py-1){
                ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                deltay=ystep_n;
            }else{
                ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                deltay=abs(ystep_p+ystep_n)/2;
            }

            JSn+= (DDmaterial[pointer_ip].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                 (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay/xstep;

            JSp+= (DDmaterial[pointer_ip].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                 (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay/xstep;
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalD2D_ISFET(double &JDn, double &JDp){

    for (int j=0; j<py; j++) {
        int pointer = (px)*(j) + (px-2);

        if(mesh[pointer].coordY>SubstrateThickness-JunctionLength && mesh[pointer].coordY<=SubstrateThickness){

            int pointer_ip = (px)*(j) + (px-1);
            int pointer_jp = (px)*(j+1) + (px-2);
            int pointer_jn = (px)*(j-1) + (px-1);
            double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
            double deltay=0;
            double ystep_p=0;
            double ystep_n=0;

            if(j==0){
                ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                deltay=ystep_p;
            }
            else if(j==py-1){
                ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                deltay=ystep_n;
            }else{
                ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                deltay=abs(ystep_p+ystep_n)/2;
            }
            JDn+= (DDmaterial[pointer].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                 (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay/xstep;

            JDp+= (DDmaterial[pointer].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                 (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay/xstep;
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalB2D_ISFET(double &JBn, double &JBp){

    for (int i=0; i<px; i++) {
        int pointer = (px)*(1) + (i);
        int pointer_jn = (px)*(0) + (i);
        int pointer_ip = (px)*(1) + (i+1);
        int pointer_in = (px)*(1) + (i-1);

        double ystep=abs(mesh[pointer].coordY-mesh[pointer_jn].coordY);

        double deltax=0;
        double xstep_p=0;
        double xstep_n=0;

        if(i==0){
            xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
            deltax=xstep_p;
        }
        else if(i==px-1){
            xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
            deltax=xstep_n;
        }else{
            xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
            xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
            deltax=abs(xstep_p+xstep_n)/2;
        }

        JBn+= (DDmaterial[pointer_jn].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_jn].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jn].phi/VT)*
             (DDmaterial[pointer_jn].u-DDmaterial[pointer].u)*deltax/ystep;

        JBp+= (DDmaterial[pointer].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_jn].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
             (DDmaterial[pointer_jn].v-DDmaterial[pointer].v)*deltax/ystep;
    }

    JBn=JBn*q0*ni_nm*VT*(-1);
    JBp=JBp*q0*ni_nm*VT;

}

void DDmodel::DD_Jcal3D(){

    switch(StructureFlag){

    case 1:
        DD_Jcal3D_PN();
        break;
    case 2:
        DD_Jcal3D_MOSFET();
        break;
    case 3:
        DD_Jcal3D_ISFET();
        break;
    case 4:
        DD_Jcal3D_1NWR();
        break;
    case 5:
        DD_Jcal3D_2NWR();
        break;
    default:
        cout << "Not appropriate StructureFlag @ DD_Jcal3D."<<endl;
        exit(0);
    }
}

double DDmodel::DD_ECSolver3D(){

    DD_loop=0;
    double errEC(0),errEC_max(0);

    do{
        DD_loop++;

        switch(StructureFlag){

        case 1:
            errEC=DD_ECTypeA3D();
            break;
        case 2:
            errEC=DD_ECTypeA3D();
            break;
        case 3:
            errEC=DD_ECTypeB3D();
            break;
        case 4:
            errEC=DD_EC1NWR3D();
            break;
        case 5:
            errEC=DD_EC2NWR3D();
            break;
        default:
            cout << "Not appropriate StructureFlag @ ECSolver3D." <<endl;
            exit(0);
        }

        if(errEC_max < errEC) {errEC_max=errEC;}

        if(DD_loop%1000==0)
        cerr <<"EC:"<< DD_loop <<"\t" <<errEC <<"\t"<<errEC_max<<endl;

    }while(errEC>SimTolEC);

    return errEC_max;

}

double DDmodel::DD_ECTypeA3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<pz-1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = DDmaterial[pointer].u;

                DDmaterial[pointer].u=DD_ECInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_ECBC3D();

    return max_val;
}

double DDmodel::DD_HCTypeA3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<pz-1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = DDmaterial[pointer].v;

                DDmaterial[pointer].v=DD_HCInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_HCBC3D();

    return max_val;
}

double DDmodel::DD_ECTypeB3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<SubstrateTOP; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = DDmaterial[pointer].u;

                DDmaterial[pointer].u=DD_ECInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_ECBC3D();

    return max_val;
}

double DDmodel::DD_HCTypeB3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<SubstrateTOP; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = DDmaterial[pointer].v;

                DDmaterial[pointer].v=DD_HCInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_HCBC3D();

    return max_val;
}

double DDmodel::DD_EC1NWR3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft1+1; i<NWRright1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom1+1; k<NWRtop1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = DDmaterial[pointer].u;

                DDmaterial[pointer].u=DD_ECInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_ECBC3D();

    return max_val;

}

double DDmodel::DD_EC2NWR3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft1+1; i<NWRright1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom1+1; k<NWRtop1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = DDmaterial[pointer].u;

                DDmaterial[pointer].u=DD_ECInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft2+1; i<NWRright2; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom2+1; k<NWRtop2; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = DDmaterial[pointer].u;

                DDmaterial[pointer].u=DD_ECInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_ECBC3D_2NWR();

    return max_val;

}

double DDmodel::DD_ECInner3D(int i, int j, int k){

    int pointer = (px)*(py)*(k) + (px)*(j) + (i);
    int pointer_ip = (px)*(py)*(k) + (px)*(j) + (i+1);
    int pointer_in = (px)*(py)*(k) + (px)*(j) + (i-1);
    int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (i);
    int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (i);
    int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (i);
    int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (i);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
    double deltaz=abs(mesh[pointer_kp].coordZ-mesh[pointer_kn].coordZ)/2;

    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
    double zstep_p=abs(mesh[pointer_kp].coordZ-mesh[pointer].coordZ);
    double zstep_n=abs(mesh[pointer_kn].coordZ-mesh[pointer].coordZ);

    double volume;

    volume=deltax*deltay*deltaz;

    double uf =DDmaterial[pointer].u;

    double Bip = (DDmaterial[pointer].mun+DDmaterial[pointer_ip].mun)/2.0*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT);
    double Bin = (DDmaterial[pointer].mun+DDmaterial[pointer_in].mun)/2.0*DD_Bern(DDmaterial[pointer_in].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_in].phi/VT);
    double Bjp = (DDmaterial[pointer].mun+DDmaterial[pointer_jp].mun)/2.0*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jp].phi/VT);
    double Bjn = (DDmaterial[pointer].mun+DDmaterial[pointer_jn].mun)/2.0*DD_Bern(DDmaterial[pointer_jn].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_jn].phi/VT);
    double Bkp = (DDmaterial[pointer].mun+DDmaterial[pointer_kp].mun)/2.0*DD_Bern(DDmaterial[pointer_kp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_kp].phi/VT);
    double Bkn = (DDmaterial[pointer].mun+DDmaterial[pointer_kn].mun)/2.0*DD_Bern(DDmaterial[pointer_kn].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_kn].phi/VT);

    double f=volume/VT*DDmaterial[pointer].r;

    double df=volume/VT*(pow(DDmaterial[pointer].v,2)*exp(-DDmaterial[pointer].phi/VT)+2*DDmaterial[pointer].v+exp(DDmaterial[pointer].phi/VT))/DDmaterial[pointer].tau/
                pow(DDmaterial[pointer].u*exp(DDmaterial[pointer].phi/VT)+DDmaterial[pointer].v*exp(-DDmaterial[pointer].phi/VT)+2,2);

    uf=((Bip*DDmaterial[pointer_ip].u/xstep_p+Bin*DDmaterial[pointer_in].u/xstep_n)*deltay*deltaz
       +(Bjp*DDmaterial[pointer_jp].u/ystep_p+Bjn*DDmaterial[pointer_jn].u/ystep_n)*deltax*deltaz
       +(Bkp*DDmaterial[pointer_kp].u/zstep_p+Bkn*DDmaterial[pointer_kn].u/zstep_n)*deltax*deltay-f+df*uf)
       /((Bip/xstep_p+Bin/xstep_n)*deltay*deltaz+(Bjp/ystep_p+Bjn/ystep_n)*deltax*deltaz+(Bkp/zstep_p+Bkn/zstep_n)*deltax*deltay+df);

    return uf;
}

double DDmodel::DD_HCSolver3D(){

    DD_loop=0;
    double errHC(0),errHC_max(0);

    do{
        DD_loop++;

        switch(StructureFlag){

        case 1:
            errHC=DD_HCTypeA3D();
            break;
        case 2:
            errHC=DD_HCTypeA3D();
            break;
        case 3:
            errHC=DD_HCTypeB3D();
            break;
        case 4:
            errHC=DD_HC1NWR3D();
            break;
        case 5:
            errHC=DD_HC2NWR3D();
            break;
        default:
            cout << "Not appropriate StructureFlag @ HCSolver3D." <<endl;
            exit(0);
        }

        if(errHC_max < errHC) {errHC_max=errHC;}

        if(DD_loop%1000==0)
        cerr <<"HC:"<< DD_loop <<"\t" <<errHC <<"\t"<<errHC_max<<endl;

    }while(errHC>SimTolHC);

    return errHC_max;
}

double DDmodel::DD_HC1NWR3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft1+1; i<NWRright1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom1+1; k<NWRtop1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = DDmaterial[pointer].v;

                DDmaterial[pointer].v=DD_HCInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_HCBC3D();

    return max_val;

}

double DDmodel::DD_HC2NWR3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft1+1; i<NWRright1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom1+1; k<NWRtop1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = DDmaterial[pointer].v;

                DDmaterial[pointer].v=DD_HCInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft2+1; i<NWRright2; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom2+1; k<NWRtop2; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = DDmaterial[pointer].v;

                DDmaterial[pointer].v=DD_HCInner3D(i,j,k);

                DDmaterial[pointer].r=DD_SRHrecomb3D(i,j,k);

                double error=VT*abs(log(DDmaterial[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    DD_HCBC3D_2NWR();

    return max_val;

}

double DDmodel::DD_HCInner3D(int i, int j, int k){

    int pointer = (px)*(py)*(k) + (px)*(j) + (i);
    int pointer_ip = (px)*(py)*(k) + (px)*(j) + (i+1);
    int pointer_in = (px)*(py)*(k) + (px)*(j) + (i-1);
    int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (i);
    int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (i);
    int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (i);
    int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (i);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
    double deltaz=abs(mesh[pointer_kp].coordZ-mesh[pointer_kn].coordZ)/2;

    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
    double zstep_p=abs(mesh[pointer_kp].coordZ-mesh[pointer].coordZ);
    double zstep_n=abs(mesh[pointer_kn].coordZ-mesh[pointer].coordZ);

    double volume;

    volume=deltax*deltay*deltaz;

    double vf =DDmaterial[pointer].v;

    double Bip = (DDmaterial[pointer].mup+DDmaterial[pointer_ip].mup)/2.0*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);
    double Bin = (DDmaterial[pointer].mup+DDmaterial[pointer_in].mup)/2.0*DD_Bern(DDmaterial[pointer_in].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);
    double Bjp = (DDmaterial[pointer].mup+DDmaterial[pointer_jp].mup)/2.0*DD_Bern(DDmaterial[pointer_jp].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);
    double Bjn = (DDmaterial[pointer].mup+DDmaterial[pointer_jn].mup)/2.0*DD_Bern(DDmaterial[pointer_jn].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);
    double Bkp = (DDmaterial[pointer].mup+DDmaterial[pointer_kp].mup)/2.0*DD_Bern(DDmaterial[pointer_kp].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);
    double Bkn = (DDmaterial[pointer].mup+DDmaterial[pointer_kn].mup)/2.0*DD_Bern(DDmaterial[pointer_kn].phi/VT-DDmaterial[pointer].phi/VT, (-1)*DDmaterial[pointer].phi/VT);

    double f=volume/VT*DDmaterial[pointer].r;

    double df=volume/VT*(pow(DDmaterial[pointer].u,2)*exp(DDmaterial[pointer].phi/VT)+2*DDmaterial[pointer].u+exp(-DDmaterial[pointer].phi/VT))/DDmaterial[pointer].tau/
                pow(DDmaterial[pointer].u*exp(DDmaterial[pointer].phi/VT)+DDmaterial[pointer].v*exp(-DDmaterial[pointer].phi/VT)+2,2);

    vf=((Bip*DDmaterial[pointer_ip].v/xstep_p+Bin*DDmaterial[pointer_in].v/xstep_n)*deltay*deltaz
       +(Bjp*DDmaterial[pointer_jp].v/ystep_p+Bjn*DDmaterial[pointer_jn].v/ystep_n)*deltax*deltaz
       +(Bkp*DDmaterial[pointer_kp].v/zstep_p+Bkn*DDmaterial[pointer_kn].v/zstep_n)*deltax*deltay-f+df*vf)
       /((Bip/xstep_p+Bin/xstep_n)*deltay*deltaz+(Bjp/ystep_p+Bjn/ystep_n)*deltax*deltaz+(Bkp/zstep_p+Bkn/zstep_n)*deltax*deltay+df);

    return vf;
}


void DDmodel::DD_Jcal3D_PN(){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0);

    DD_JcalS3D_PN(Current_Sn,Current_Sp);
    DD_JcalD3D_PN(Current_Dn,Current_Dp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific<<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp<<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<endl;
    output1.close();
}

void DDmodel::DD_JcalS3D_PN(double &JSn,double &JSp){

    for (int j=0; j<py; j++) {
        for (int k=1; k<pz-1; k++) {
            int pointer = (px)*(py)*(k) + (px)*(j) + (0);
            int pointer_ip= (px)*(py)*(k) + (px)*(j) + (1);
            int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (0);
            int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (0);
            int pointer_kp = (px)*(py)*(k+1) + (px)*(0) + (0);
            int pointer_kn = (px)*(py)*(k-1) + (px)*(0) + (0);

            double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
            double deltay=0;
            double deltaz=0;

            if(j==0){
                deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
            }
            else if(j==py-1){
                deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
            }else {
                deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
            }

            if(k==0){
                deltaz=abs(mesh[pointer_kp].coordZ-mesh[pointer].coordZ);
            }
            else if(k==pz-1){
                deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer].coordZ);
            }else {
                deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;
            }

            JSn+= (DDmaterial[pointer].mun+DDmaterial[pointer_ip].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                 (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay*deltaz/xstep;

            JSp+= (DDmaterial[pointer].mup+DDmaterial[pointer_ip].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer_ip].phi/VT)*
                 (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay*deltaz/xstep;
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalD3D_PN(double &JDn, double &JDp){

    for (int j=0; j<py; j++) {
        for (int k=0; k<pz; k++) {
            int pointer = (px)*(py)*(k) + (px)*(j) + (px-2);
            int pointer_ip= (px)*(py)*(k) + (px)*(j) + (px-1);
            int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (px-1);
            int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (px-1);
            int pointer_kp = (px)*(py)*(k+1) + (px)*(0) + (px-1);
            int pointer_kn = (px)*(py)*(k-1) + (px)*(0) + (px-1);

            double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
            double deltay=0;
            double deltaz=0;

            if(j==0){
                deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
            }
            else if(j==py-1){
                deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
            }else {
                deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
            }

            if(k==0){
                deltaz=abs(mesh[pointer_kp].coordZ-mesh[pointer].coordZ);
            }
            else if(k==pz-1){
                deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer].coordZ);
            }else {
                deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;
            }

            JDn+= (DDmaterial[pointer].mun+DDmaterial[pointer_ip].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                 (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay*deltaz/xstep;

            JDp+= (DDmaterial[pointer].mup+DDmaterial[pointer_ip].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer_ip].phi/VT)*
                 (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay*deltaz/xstep;
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}


void DDmodel::DD_Jcal3D_MOSFET(){

    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    DD_JcalS3D_MOSFET(Current_Sn,Current_Sp);
    DD_JcalD3D_MOSFET(Current_Dn,Current_Dp);
    DD_JcalB3D_MOSFET(Current_Bn,Current_Bp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific
          <<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
          <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"
          <<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;
    output1.close();
}

void DDmodel::DD_JcalS3D_MOSFET(double &JSn, double &JSp){

    for (int i=0; i<px; i++) {
        for (int j=0; j<py; j++) {

            int pointer = (px)*(py)*(0) + (px)*(j) + (i);

            if(mesh[pointer].coordX<=JunctionLength){
                int pointer_kp = (px)*(py)*(1) + (px)*(j) + (i);
                int pointer_in = (px)*(py)*(0) + (px)*(j) + (i-1);
                int pointer_ip = (px)*(py)*(0) + (px)*(j) + (i+1);
                int pointer_jp = (px)*(py)*(0) + (px)*(j+1) + (i);
                int pointer_jn = (px)*(py)*(0) + (px)*(j-1) + (i);

                double zstep=abs(mesh[pointer].coordZ-mesh[pointer_kp].coordZ);
                double deltax=0;
                double deltay=0;

                if(i==0){
                    deltax=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
                }
                else if(i==px-1){
                    deltax=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
                }else{
                    deltax=abs(mesh[pointer_in].coordX-mesh[pointer_ip].coordX)/2;
                }

                if(j==0){
                    deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                }
                else if(j==py-1){
                    deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                }else{
                    deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                }

                JSn+= (DDmaterial[pointer_kp].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_kp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_kp].phi/VT)*
                     (DDmaterial[pointer_kp].u-DDmaterial[pointer].u)*deltax*deltay/zstep;

                JSp+= (DDmaterial[pointer_kp].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_kp].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                     (DDmaterial[pointer_kp].v-DDmaterial[pointer].v)*deltax*deltay/zstep;

            }
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalD3D_MOSFET(double &JDn, double &JDp){

    for (int i=0; i<px; i++) {
        for (int j=0; j<py; j++) {

            int pointer = (px)*(py)*(0) + (px)*(j) + (i);

            if(mesh[pointer].coordX>=lx-JunctionLength){
                int pointer_kp = (px)*(py)*(1) + (px)*(j) + (i);
                int pointer_in = (px)*(py)*(0) + (px)*(j) + (i-1);
                int pointer_ip = (px)*(py)*(0) + (px)*(j) + (i+1);
                int pointer_jp = (px)*(py)*(0) + (px)*(j+1) + (i);
                int pointer_jn = (px)*(py)*(0) + (px)*(j-1) + (i);

                double zstep=abs(mesh[pointer].coordZ-mesh[pointer_kp].coordZ);
                double deltax=0;
                double deltay=0;

                if(i==0){
                    deltax=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
                }
                else if(i==px-1){
                    deltax=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
                }else{
                    deltax=abs(mesh[pointer_in].coordX-mesh[pointer_ip].coordX)/2;
                }

                if(j==0){
                    deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                }
                else if(j==py-1){
                    deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                }else{
                    deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                }

                JDn+= (DDmaterial[pointer_kp].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_kp].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_kp].phi/VT)*
                     (DDmaterial[pointer_kp].u-DDmaterial[pointer].u)*deltax*deltay/zstep;

                JDp+= (DDmaterial[pointer_kp].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_kp].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                     (DDmaterial[pointer_kp].v-DDmaterial[pointer].v)*deltax*deltay/zstep;
            }
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalB3D_MOSFET(double &JBn, double &JBp){

    for (int i=0; i<px; i++) {
        for (int j=0; j<py; j++) {

            int pointer = (px)*(py)*(pz-1) + (px)*(j) + (i);
            int pointer_kn = (px)*(py)*(px-2) + (px)*(j) + (i);
            int pointer_in = (px)*(py)*(pz-1) + (px)*(j) + (i-1);
            int pointer_ip = (px)*(py)*(pz-1) + (px)*(j) + (i+1);
            int pointer_jp = (px)*(py)*(pz-1) + (px)*(j+1) + (i);
            int pointer_jn = (px)*(py)*(pz-1) + (px)*(j-1) + (i);

            double zstep=abs(mesh[pointer].coordZ-mesh[pointer_kn].coordZ);
            double deltax=0;
            double deltay=0;

            if(i==0){
                deltax=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
            }
            else if(i==px-1){
                deltax=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
            }else{
                deltax=abs(mesh[pointer_in].coordX-mesh[pointer_ip].coordX)/2;
            }

            if(j==0){
                deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
            }
            else if(j==py-1){
                deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
            }else{
                deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
            }

            JBn+= (DDmaterial[pointer_kn].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_kn].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_kn].phi/VT)*
                 (DDmaterial[pointer_kn].u-DDmaterial[pointer].u)*deltax*deltay/zstep;

            JBp+= (DDmaterial[pointer_kn].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_kn].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                 (DDmaterial[pointer_kn].v-DDmaterial[pointer].v)*deltax*deltay/zstep;
        }
    }

    JBn=JBn*q0*ni_nm*VT*(-1);
    JBp=JBp*q0*ni_nm*VT;
}

void DDmodel::DD_Jcal3D_ISFET(){

    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    DD_JcalS3D_ISFET(Current_Sn,Current_Sp);
    DD_JcalD3D_ISFET(Current_Dn,Current_Dp);
    DD_JcalB3D_ISFET(Current_Bn,Current_Bp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific
          <<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
          <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"
          <<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;
    output1.close();
}

void DDmodel::DD_JcalS3D_ISFET(double &JSn, double &JSp){

    for (int j=0; j<py; j++) {
        for (int k=0; k<pz; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (0);

            if(mesh[pointer].coordZ>SubstrateThickness-JunctionDepth && mesh[pointer].coordZ<=SubstrateThickness ){

                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (0);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (0);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (0);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (0);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=0;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;;

                if(j==0){
                    deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                }
                else if(j==py-1){
                    deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                }else{
                    deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                }

                JSn+= (DDmaterial[pointer_ip].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                     (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay*deltaz/xstep;

                JSp+= (DDmaterial[pointer_ip].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                     (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalD3D_ISFET(double &JDn, double &JDp){

    for (int j=0; j<py; j++) {
        for (int k=0; k<pz; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (px-2);

            if(mesh[pointer].coordZ>SubstrateThickness-JunctionDepth && mesh[pointer].coordZ<=SubstrateThickness ){

                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (px-1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (px-2);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (px-2);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (px-2);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (px-2);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=0;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;;

                if(j==0){
                    deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                }
                else if(j==py-1){
                    deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                }else{
                    deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                }

                JDn+= (DDmaterial[pointer_ip].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                     (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltaz*deltay/xstep;

                JDp+= (DDmaterial[pointer_ip].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                     (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltaz*deltay/xstep;
            }
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalB3D_ISFET(double &JBn, double &JBp){

    for (int i=0; i<px; i++) {
        for (int j=0; j<py; j++) {

            int pointer = (px)*(py)*(pz-1) + (px)*(j) + (i);
            int pointer_kn = (px)*(py)*(px-2) + (px)*(j) + (i);
            int pointer_in = (px)*(py)*(pz-1) + (px)*(j) + (i-1);
            int pointer_ip = (px)*(py)*(pz-1) + (px)*(j) + (i+1);
            int pointer_jp = (px)*(py)*(pz-1) + (px)*(j+1) + (i);
            int pointer_jn = (px)*(py)*(pz-1) + (px)*(j-1) + (i);

            double zstep=abs(mesh[pointer].coordZ-mesh[pointer_kn].coordZ);
            double deltax=0;
            double deltay=0;

            if(i==0){
                deltax=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
            }
            else if(i==px-1){
                deltax=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
            }else{
                deltax=abs(mesh[pointer_in].coordX-mesh[pointer_ip].coordX)/2;
            }

            if(j==0){
                deltay=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
            }
            else if(j==py-1){
                deltay=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
            }else{
                deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
            }

            JBn+= (DDmaterial[pointer_kn].mun+DDmaterial[pointer].mun)/2*DD_Bern(DDmaterial[pointer_kn].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_kn].phi/VT)*
                 (DDmaterial[pointer_kn].u-DDmaterial[pointer].u)*deltax*deltay/zstep;

            JBp+= (DDmaterial[pointer_kn].mup+DDmaterial[pointer].mup)/2*DD_Bern(DDmaterial[pointer_kn].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                 (DDmaterial[pointer_kn].v-DDmaterial[pointer].v)*deltax*deltay/zstep;
        }
    }

    JBn=JBn*q0*ni_nm*VT*(-1);
    JBp=JBp*q0*ni_nm*VT;
}


void DDmodel::DD_Jcal3D_1NWR(){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    DD_JcalS3D_NWR1(Current_Sn,Current_Sp);
    DD_JcalD3D_NWR1(Current_Dn,Current_Dp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific<<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
           <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"<<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;

    output1.close();
}

void DDmodel::DD_Jcal3D_2NWR(){
    ofstream  output1;
    double Current_Sn1(0),Current_Sp1(0),Current_Dn1(0),Current_Dp1(0);
    double Current_Sn2(0),Current_Sp2(0),Current_Dn2(0),Current_Dp2(0);

    DD_JcalS3D_NWR1(Current_Sn1,Current_Sp1);
    DD_JcalD3D_NWR1(Current_Dn1,Current_Dp1);

    DD_JcalS3D_NWR2(Current_Sn2,Current_Sp2);
    DD_JcalD3D_NWR2(Current_Dn2,Current_Dp2);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<(NWRCentery1+NWRCentery2)/2+ShiftDistanceY<<"\t"<<scientific
          <<Current_Sn1<<"\t"<<Current_Sp1<<"\t"<<Current_Dn1<<"\t"<<Current_Dp1<<"\t"<<Current_Sn1+Current_Sp1<<"\t"<<Current_Dn1+Current_Dp1<<"\t"
          <<Current_Sn2<<"\t"<<Current_Sp2<<"\t"<<Current_Dn2<<"\t"<<Current_Dp2<<"\t"<<Current_Sn2+Current_Sp2<<"\t"<<Current_Dn2+Current_Dp2<<"\t"<<endl;

    output1.close();
}

void DDmodel::DD_JcalS3D_NWR1(double &JSn, double &JSp){

    for (int j=NWRleft1; j<NWRright1; j++) {
        for (int k=NWRbottom1; k<NWRtop1; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (0);

            if(DDmaterial[pointer].Type==1){
                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (1);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (1);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (1);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (1);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;

                JSn+= (DDmaterial[pointer].mun+DDmaterial[pointer_ip].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                     (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay*deltaz/xstep;

                JSp+=(DDmaterial[pointer].mup+DDmaterial[pointer_ip].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                    (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalD3D_NWR1(double &JDn, double &JDp){

    for (int j=NWRleft1; j<NWRright1; j++) {
        for (int k=NWRbottom1; k<NWRtop1; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (px-2);
            if(DDmaterial[pointer].Type==1){
                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (px-1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (px-1);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (px-1);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (px-1);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (px-1);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;

                JDn+= (DDmaterial[pointer].mun+DDmaterial[pointer_ip].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                     (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay*deltaz/xstep;

                JDp+=(DDmaterial[pointer].mup+DDmaterial[pointer_ip].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                    (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalS3D_NWR2(double &JSn, double &JSp){

    for (int j=NWRleft1; j<NWRright1; j++) {
        for (int k=NWRbottom1; k<NWRtop1; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (0);

            if(DDmaterial[pointer].Type==1){
                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (1);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (1);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (1);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (1);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;

                JSn+= (DDmaterial[pointer].mun+DDmaterial[pointer_ip].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                     (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay*deltaz/xstep;

                JSp+=(DDmaterial[pointer].mup+DDmaterial[pointer_ip].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                    (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::DD_JcalD3D_NWR2(double &JDn, double &JDp){

    for (int j=NWRleft1; j<NWRright1; j++) {
        for (int k=NWRbottom1; k<NWRtop1; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (px-2);
            if(DDmaterial[pointer].Type==1){
                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (px-1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (px-1);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (px-1);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (px-1);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (px-1);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;

                JDn+= (DDmaterial[pointer].mun+DDmaterial[pointer_ip].mun)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, DDmaterial[pointer_ip].phi/VT)*
                     (DDmaterial[pointer_ip].u-DDmaterial[pointer].u)*deltay*deltaz/xstep;

                JDp+=(DDmaterial[pointer].mup+DDmaterial[pointer_ip].mup)/2*DD_Bern(DDmaterial[pointer_ip].phi/VT-DDmaterial[pointer].phi/VT, -DDmaterial[pointer].phi/VT)*
                    (DDmaterial[pointer_ip].v-DDmaterial[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::DD_FindNWRBC1(){

    //find center point
    for(int j=0;j<py;j++){
        int pointer = (py)*(px)*(pz/2) + (px)*(j) + (px/2);
        if(abs(mesh[pointer].coordY-NWRCentery1)<0.01)
            NWRcenteryP1=j;
    }

    for(int k=0;k<pz;k++){
        int pointer = (py)*(px)*(k) + (px)*(py/2) + (px/2);
        if(abs(mesh[pointer].coordZ-NWRCenterz1)<0.01)
            NWRcenterzP1=k;
    }

    //cerr << "(NWRcenteryP1,NWRcenterzP1)="<<"("<<NWRcenteryP1<<","<<NWRcenterzP1<<")"<<endl;

    //find boundary for NWR:tox
    for(int k=NWRcenterzP1;k<pz;k++){

        int pointer = (py)*(px)*(k) + (px)*(NWRcenteryP1) + (px/2);
        int pointer_kp = (py)*(px)*(k+1) + (px)*(NWRcenteryP1) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_kp].Type==2){
            NWRtop1=k;
            break;
        }
    }

    for(int k=NWRcenterzP1;k>0;k--){

        int pointer = (py)*(px)*(k) + (px)*(NWRcenteryP1) + (px/2);
        int pointer_kn = (py)*(px)*(k-1) + (px)*(NWRcenteryP1) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_kn].Type==2){
            NWRbottom1=k;
            break;
        }
    }

    //cerr << "(NWRbottom1,NWRtop1)="<<"("<<NWRbottom1<<","<<NWRtop1<<")"<<endl;

    for(int j=NWRcenteryP1;j<py;j++){

        int pointer = (py)*(px)*(NWRcenterzP1) + (px)*(j) + (px/2);
        int pointer_jp = (py)*(px)*(NWRcenterzP1) + (px)*(j+1) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_jp].Type==2){
            NWRright1=j;
            break;
        }
    }

    for(int j=NWRcenteryP1;j<py;j--){

        int pointer = (py)*(px)*(NWRcenterzP1) + (px)*(j) + (px/2);
        int pointer_jn = (py)*(px)*(NWRcenterzP1) + (px)*(j-1) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_jn].Type==2){
            NWRleft1=j;
            break;
        }
    }
    //cerr << "(NWRleft1,NWRright1)="<<"("<<NWRleft1<<","<<NWRright1<<")"<<endl;

}

void DDmodel::DD_FindNWRBC2(){

    //find center point
    for(int j=0;j<py;j++){
        int pointer = (py)*(px)*(pz/2) + (px)*(j) + (px/2);
        if(abs(mesh[pointer].coordY-NWRCentery2)<0.01)
            NWRcenteryP2=j;
    }

    for(int k=0;k<pz;k++){
        int pointer = (py)*(px)*(k) + (px)*(py/2) + (px/2);
        if(abs(mesh[pointer].coordZ-NWRCenterz2)<0.01)
            NWRcenterzP2=k;
    }

    //cerr << "(NWRcenteryP2,NWRcenterzP2)="<<"("<<NWRcenteryP2<<","<<NWRcenterzP2<<")"<<endl;

    //find boundary for NWR:tox
    for(int k=NWRcenterzP2;k<pz;k++){

        int pointer = (py)*(px)*(k) + (px)*(NWRcenteryP2) + (px/2);
        int pointer_kp = (py)*(px)*(k+1) + (px)*(NWRcenteryP2) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_kp].Type==2){
            NWRtop2=k;
            break;
        }
    }

    for(int k=NWRcenterzP2;k>0;k--){

        int pointer = (py)*(px)*(k) + (px)*(NWRcenteryP2) + (px/2);
        int pointer_kn = (py)*(px)*(k-1) + (px)*(NWRcenteryP2) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_kn].Type==2){
            NWRbottom2=k;
            break;
        }
    }
    //cerr << "(NWRbottom2,NWRtop2)="<<"("<<NWRbottom2<<","<<NWRtop2<<")"<<endl;


    for(int j=NWRcenteryP2;j<py;j++){

        int pointer = (py)*(px)*(NWRcenterzP2) + (px)*(j) + (px/2);
        int pointer_jp = (py)*(px)*(NWRcenterzP2) + (px)*(j+1) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_jp].Type==2){
            NWRright2=j;
            break;
        }
    }

    for(int j=NWRcenteryP2;j<py;j--){

        int pointer = (py)*(px)*(NWRcenterzP2) + (px)*(j) + (px/2);
        int pointer_jn = (py)*(px)*(NWRcenterzP2) + (px)*(j-1) + (px/2);

        if(DDmaterial[pointer].Type==1 && DDmaterial[pointer_jn].Type==2){
            NWRleft2=j;
            break;
        }
    }
    //cerr << "(NWRleft2,NWRright2)="<<"("<<NWRleft2<<","<<NWRright2<<")"<<endl;
}

void DDmodel::DD_AddDot3D(double DotXCenter, double DotYCenter, int Flag){

    double volume=0;
    double DotZCenter=SubstrateThickness+BOX+2*NWRradiusz+Tox+AnalyteRadius+ReceptorLength;

    for(int k=0;k<pz;k++){
        for(int j=0;j<py;j++){
            for(int i=0;i<px;i++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                if(pow(mesh[pointer].coordX-DotXCenter,2)+pow(mesh[pointer].coordY-DotYCenter,2)+pow(mesh[pointer].coordZ-DotZCenter,2)<=AnalyteRadius*AnalyteRadius){
                    int pointer_ip = (py)*(px)*(k) + (px)*(j) + (i+1);
                    int pointer_in = (py)*(px)*(k) + (px)*(j) + (i-1);
                    int pointer_jp = (py)*(px)*(k) + (px)*(j+1) + (i);
                    int pointer_jn = (py)*(px)*(k) + (px)*(j-1) + (i);
                    int pointer_kp = (py)*(px)*(k+1) + (px)*(j) + (i);
                    int pointer_kn = (py)*(px)*(k-1) + (px)*(j) + (i);
                    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
                    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
                    double deltaz=abs(mesh[pointer_kp].coordZ-mesh[pointer_kn].coordZ)/2;

                    volume=volume+deltax*deltay*deltaz;

                    DDmaterial[pointer].Type=Flag;
                    DDmaterial[pointer].k=AnalytePermittivity;
                }
            }
        }
    }

    for(int k=0;k<pz;k++){
        for(int j=0;j<py;j++){
            for(int i=0;i<px;i++){

                int pointer= (py)*(px)*(k) + (px)*(j) + (i);

                if(pow(mesh[pointer].coordX-DotXCenter,2)+pow(mesh[pointer].coordY-DotYCenter,2)+pow(mesh[pointer].coordZ-DotZCenter,2)<=AnalyteRadius*AnalyteRadius){
                    DDmaterial[pointer].rho=AnalyteValence*q0/volume; //q0/volumez
                    //rho = [C/nm3]
                }
            }
        }
    }
}


void DDmodel::DD_AddDotString3D(double Yshift){

    double DotDistance=NWRLength/DotNumber;

    for(int i=0;i<DotNumber;i++){
        DD_AddDot3D(JunctionLength+DotDistance/2+DotDistance*i,(NWRCentery1+NWRCentery2)/2-Yshift,i+1000);
    }
}

