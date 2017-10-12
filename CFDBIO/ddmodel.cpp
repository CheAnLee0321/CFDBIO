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
void DDmodel::DDmodelParameterSet(){

    BernoulliX();

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

void DDmodel::DDInitialGuess2D(){

    material=new Semiconductor [L];
    DDInitialize();

    switch(StructureFlag){

    case 1:
        DDInitialGuessPNJunction2D();
        break;
    case 2:
        DDInitialGuessMOSFET2D();
        break;
    case 3:
        DDInitialGuessISFET2D();
        break;
    default:
        cout << "Undifined Device Structure @ DDInitialGuess2D." << endl;
        exit(0);
    }

}

void DDmodel::DDInitialGuessPNJunction2D(){

#pragma omp parallel for
    for(int i=0;i<L;i++){
       material[i].k=Si_permi;
       material[i].Type=1;
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer = (px)*(j) + (i);

            //setup P+ Drain
            //P+
            material[pointer].dop=-Nai;
            material[pointer].phi=(volD-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
            material[pointer].u=exp((-1)*volD/VT);
            material[pointer].v=exp(volD/VT);
            material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
            material[pointer].mup=mupCal(Tamb, 0, 1);
            material[pointer].tau=tauPCal(0);
            material[pointer].r=SRHrecomb2D(i,j);

            //setup N+ Source
            if(mesh[pointer].coordX<lx/2){
                //N+
                material[pointer].dop=Ndi;
                material[pointer].phi=(volS+VT*log(0.5*Ndi+sqrt(pow(0.5*Ndi,2)+1)));
                material[pointer].u=exp((-1)*volS/VT);
                material[pointer].v=exp(volS/VT);
                material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                material[pointer].mup=mupCal(Tamb, 0, 1);
                material[pointer].tau=tauPCal(0);
                material[pointer].r=SRHrecomb2D(i,j);
            }
        }
    }
}

void DDmodel::DDInitialGuessMOSFET2D(){

#pragma omp parallel for
    for(int i=0;i<L;i++){
       material[i].k=Si_permi;
       material[i].Type=1;
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer = (px)*(j) + (i);

            //setup P+
            //P+
            material[pointer].dop=-Nai;
            material[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
            material[pointer].u=exp((-1)*volB/VT);
            material[pointer].v=exp(volB/VT);
            material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
            material[pointer].mup=mupCal(Tamb, 0, 1);
            material[pointer].tau=tauPCal(0);
            material[pointer].r=SRHrecomb2D(i,j);

            //setup N+ Source
            if(mesh[pointer].coordX<=JunctionLength){
                if(mesh[pointer].coordY<=JunctionDepth){
                    //N+
                    material[pointer].dop=NdPlusi;
                    material[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                    material[pointer].u=exp((-1)*volS/VT);
                    material[pointer].v=exp(volS/VT);
                    material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, 1);
                    material[pointer].tau=tauPCal(0);
                    material[pointer].r=SRHrecomb2D(i,j);
                }
            }

            //setup N+ Drain
            if(mesh[pointer].coordX>=lx-JunctionLength){
                if(mesh[pointer].coordY<=JunctionDepth){
                    //N+
                    material[pointer].dop=NdPlusi;
                    material[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                    material[pointer].u=exp((-1)*volD/VT);
                    material[pointer].v=exp(volD/VT);
                    material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, 1);
                    material[pointer].tau=tauPCal(0);
                    material[pointer].r=SRHrecomb2D(i,j);
                }
            }
        }
    }
}

void DDmodel::DDInitialGuessISFET2D(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer = (px)*(j) + (i);

            //substrate
            if(mesh[pointer].coordY<=SubstrateThickness ){
                //P+
                material[pointer].Type=1;
                material[pointer].k=Si_permi;
                material[pointer].dop=-Nai;
                material[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                material[pointer].u=exp((-1)*volB/VT);
                material[pointer].v=exp(volB/VT);
                material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                material[pointer].mup=mupCal(Tamb, 0, 1);
                material[pointer].tau=tauPCal(0);
                material[pointer].r=SRHrecomb2D(i,j);
            }

            //SD
            if(mesh[pointer].coordY>SubstrateThickness-JunctionDepth && mesh[pointer].coordY<=SubstrateThickness ){

                if(mesh[pointer].coordX<JunctionLength){
                    //N+
                    material[pointer].Type=1;
                    material[pointer].k=Si_permi;
                    material[pointer].dop=NdPlusi;
                    material[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                    material[pointer].u=exp((-1)*volS/VT);
                    material[pointer].v=exp(volS/VT);
                    material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, 1);
                    material[pointer].tau=tauPCal(0);
                    material[pointer].r=SRHrecomb2D(i,j);

                }
                if(mesh[pointer].coordX>lx-JunctionLength){
                    //N+
                    material[pointer].Type=1;
                    material[pointer].k=Si_permi;
                    material[pointer].dop=NdPlusi;
                    material[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                    material[pointer].u=exp((-1)*volD/VT);
                    material[pointer].v=exp(volD/VT);
                    material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, 1);
                    material[pointer].tau=tauPCal(0);
                    material[pointer].r=SRHrecomb2D(i,j);
                }
            }

            //oxide
            if(mesh[pointer].coordY>SubstrateThickness && mesh[pointer].coordY<=SubstrateThickness+Tox ){
                material[pointer].Type=2;
                material[pointer].k=SiO2_permi;
                material[pointer].phi=volG;
                material[pointer].dop=0;
                material[pointer].u=exp((-1)*volG/VT);
                material[pointer].v=exp(volG/VT);
                material[pointer].mun=0;
                material[pointer].mup=0;
                material[pointer].tau=0;
                material[pointer].r=0;
            }

            //electrolyte
            if(mesh[pointer].coordY>SubstrateThickness+Tox){
                material[pointer].Type=3;
                material[pointer].k=Water_permi;
                material[pointer].phi=volG;
                material[pointer].dop=0;
                material[pointer].u=exp((-1)*volG/VT);
                material[pointer].v=exp(volG/VT);
                material[pointer].mun=0.12*1e18;
                material[pointer].mup=0.12*1e18;
                material[pointer].tau=tauNCal(Si_ni_cm);
                material[pointer].r=0;
            }

        }
    }

    //find boundary for sub:tox, tox:eletrolyte
    for(int j=0;j<py-1;j++){

        int pointer = (px)*(j) + (px/2);
        int pointer_jn = (px)*(j-1) + (px/2);
        int pointer_jp = (px)*(j+1) + (px/2);

        if(material[pointer].Type==1 && material[pointer_jp].Type==2){
            SubstrateTOP=j;
        }

        if(material[pointer].Type==3 && material[pointer_jn].Type==2){
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

void DDmodel::DDInitialGuess3D(){

    material=new Semiconductor [L];
    DDInitialize();

    switch(StructureFlag){

    case 1:
        DDInitialGuessPNJunction3D();
        break;
    case 2:
        DDInitialGuessMOSFET3D();
        break;
    case 3:
        DDInitialGuessISFET3D();
        break;
    case 4:
        DDInitialGuess1NWR3D();
        FindNWRBC1();
        break;
    case 5:
        DDInitialGuess2NWR3D();
        FindNWRBC1();
        FindNWRBC2();
        break;
    default:
        cout << "Undifined Device Structure @ DDInitialGuess3D." << endl;
        exit(0);
    }

}

void DDmodel::DDInitialGuessPNJunction3D(){

    StructureFlag=1;

#pragma omp parallel for
    for(int i=0;i<L;i++){
       material[i].k=Si_permi;
       material[i].Type=1;
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //setup P+ Drain
                //P+
                material[pointer].dop=-Nai;
                material[pointer].phi=(volD-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                material[pointer].u=exp((-1)*volD/VT);
                material[pointer].v=exp(volD/VT);
                material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                material[pointer].mup=mupCal(Tamb, 0, 1);
                material[pointer].tau=tauPCal(0);
                material[pointer].r=SRHrecomb2D(i,j);

                //setup N+ Source
                if(mesh[pointer].coordX<lx/2){
                    //N+
                    material[pointer].dop=Ndi;
                    material[pointer].phi=(volS+VT*log(0.5*Ndi+sqrt(pow(0.5*Ndi,2)+1)));
                    material[pointer].u=exp((-1)*volS/VT);
                    material[pointer].v=exp(volS/VT);
                    material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, 1);
                    material[pointer].tau=tauPCal(0);
                    material[pointer].r=SRHrecomb2D(i,j);
                }
            }
        }
    }
}

void DDmodel::DDInitialGuessMOSFET3D(){

    //StructureFlag=2;

#pragma omp parallel for
    for(int i=0;i<L;i++){
       material[i].k=Si_permi;
       material[i].Type=1;
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //setup P+
                //P+
                material[pointer].dop=-Nai;
                material[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                material[pointer].u=exp((-1)*volB/VT);
                material[pointer].v=exp(volB/VT);
                material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                material[pointer].mup=mupCal(Tamb, 0, 1);
                material[pointer].tau=tauPCal(0);
                material[pointer].r=SRHrecomb2D(i,j);

                //setup N+ Source
                if(mesh[pointer].coordX<=JunctionLength){
                    if(mesh[pointer].coordZ<=JunctionDepth){
                        //N+
                        material[pointer].dop=NdPlusi;
                        material[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                        material[pointer].u=exp((-1)*volS/VT);
                        material[pointer].v=exp(volS/VT);
                        material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                        material[pointer].mup=mupCal(Tamb, 0, 1);
                        material[pointer].tau=tauPCal(0);
                        material[pointer].r=SRHrecomb2D(i,j);
                    }
                }

                //setup N+ Drain
                if(mesh[pointer].coordX>=lx-JunctionLength){
                    if(mesh[pointer].coordZ<=JunctionDepth){
                        //N+
                        material[pointer].dop=NdPlusi;
                        material[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                        material[pointer].u=exp((-1)*volD/VT);
                        material[pointer].v=exp(volD/VT);
                        material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                        material[pointer].mup=mupCal(Tamb, 0, 1);
                        material[pointer].tau=tauPCal(0);
                        material[pointer].r=SRHrecomb2D(i,j);
                    }
                }
            }
        }
    }
}

void DDmodel::DDInitialGuessISFET3D(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //substrate
                if(mesh[pointer].coordZ<=SubstrateThickness ){
                    //P+
                    material[pointer].Type=1;
                    material[pointer].k=Si_permi;
                    material[pointer].dop=-Nai;
                    material[pointer].phi=(volB-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                    material[pointer].u=exp((-1)*volB/VT);
                    material[pointer].v=exp(volB/VT);
                    material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, 1);
                    material[pointer].tau=tauPCal(0);
                    material[pointer].r=SRHrecomb2D(i,j);
                }

                //SD
                if(mesh[pointer].coordZ>SubstrateThickness-JunctionDepth && mesh[pointer].coordZ<=SubstrateThickness ){

                    if(mesh[pointer].coordX<JunctionLength){
                        //N+
                        material[pointer].Type=1;
                        material[pointer].k=Si_permi;
                        material[pointer].dop=NdPlusi;
                        material[pointer].phi=(volS+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                        material[pointer].u=exp((-1)*volS/VT);
                        material[pointer].v=exp(volS/VT);
                        material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                        material[pointer].mup=mupCal(Tamb, 0, 1);
                        material[pointer].tau=tauPCal(0);
                        material[pointer].r=SRHrecomb2D(i,j);
                    }
                    if(mesh[pointer].coordX>lx-JunctionLength){
                        //N+
                        material[pointer].Type=1;
                        material[pointer].k=Si_permi;
                        material[pointer].dop=NdPlusi;
                        material[pointer].phi=(volD+VT*log(0.5*NdPlusi+sqrt(pow(0.5*NdPlusi,2)+1)));
                        material[pointer].u=exp((-1)*volD/VT);
                        material[pointer].v=exp(volD/VT);
                        material[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                        material[pointer].mup=mupCal(Tamb, 0, 1);
                        material[pointer].tau=tauPCal(0);
                        material[pointer].r=SRHrecomb2D(i,j);
                    }
                }

                //oxide
                if(mesh[pointer].coordZ>SubstrateThickness && mesh[pointer].coordZ<=SubstrateThickness+Tox ){
                    material[pointer].Type=2;
                    material[pointer].k=SiO2_permi;
                    material[pointer].phi=volG;
                    material[pointer].dop=0;
                    material[pointer].u=exp((-1)*volG/VT);
                    material[pointer].v=exp(volG/VT);
                    material[pointer].mun=0;
                    material[pointer].mup=0;
                    material[pointer].tau=0;
                    material[pointer].r=0;
                }

                //electrolyte
                if(mesh[pointer].coordZ>SubstrateThickness+Tox){
                    material[pointer].Type=3;
                    material[pointer].k=Water_permi;
                    material[pointer].phi=volG;
                    material[pointer].dop=0;
                    material[pointer].u=exp((-1)*volG/VT);
                    material[pointer].v=exp(volG/VT);
                    material[pointer].mun=0.12*1e18;
                    material[pointer].mup=0.12*1e18;
                    material[pointer].tau=tauNCal(Si_ni_cm);
                    material[pointer].r=0;
                }
            }
        }
    }

    //find boundary for sub:tox, tox:eletrolyte
    for(int k=0;k<pz-1;k++){

        int pointer = (py)*(px)*(k) + (px)*(py/2) + (px/2);
        int pointer_kn = (py)*(px)*(k-1) + (px)*(py/2) + (px/2);
        int pointer_kp = (py)*(px)*(k+1) + (px)*(py/2) + (px/2);

        if(material[pointer].Type==1 && material[pointer_kp].Type==2){
            SubstrateTOP=k;
            //cerr << material[pointer_jp].Type <<'\t'<<j<<endl;
        }

        if(material[pointer].Type==3 && material[pointer_kn].Type==2){
            ElectrolyteBottom=k;
            //cerr << material[pointer_jn].Type <<'\t'<<j<<endl;
        }
    }

    if(SubstrateTOP==0||ElectrolyteBottom==0){
        cerr << "Boundary undifined."<<endl;
        cerr << "SubstrateTOP=" <<SubstrateTOP <<endl;
        cerr << "ElectrolyteBottom="<<ElectrolyteBottom<<endl;
    }
}

void DDmodel::DDInitialGuess1NWR3D(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //substrate
                if(mesh[pointer].coordZ<=SubstrateThickness ){
                    material[pointer].Type=4;
                    material[pointer].k=Si_permi;
                    material[pointer].dop=0;
                    material[pointer].phi=volB;
                    material[pointer].u=exp((-1)*volB/VT);
                    material[pointer].v=exp(volB/VT);
                    material[pointer].mun=munCal(Tamb, 0, material[pointer].Type); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, material[pointer].Type);
                    material[pointer].tau=tauPCal(0);
                    material[pointer].r=SRHrecomb3D(i,j,k);
                    material[pointer].rho=0;
                }

                //electrolyte
                if(mesh[pointer].coordZ>SubstrateThickness+BOX){
                    material[pointer].Type=3;
                    material[pointer].k=Water_permi;
                    //n+
                    //material[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2-VT*log(0.5*Ndi_SD+sqrt(pow(0.5*Ndi_SD,2)+1)))+volG;
                    //P+
                    //material[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2+VT*log(0.5*Nai_SD+sqrt(pow(0.5*Nai_SD,2)+1)))+volG;
                    //metal
                    //phi_new[pointer]=(Si_chi+(Si_Eg)/2)-(wfG)+volG;
                    //Intrinsic
                    material[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2)+volG;
                    material[pointer].dop=0;
                    material[pointer].u=exp((-1)*volG/VT);
                    material[pointer].v=exp(volG/VT);
                    material[pointer].mun=munCal(Tamb, 0, material[pointer].Type); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, material[pointer].Type);
                    material[pointer].tau=tauNCal(Si_ni_cm);
                    material[pointer].r=0;
                    material[pointer].rho=0;
                }

                //box oxide
                //Protruding
                //if(mesh[pointer].coordZ>SubstrateThickness && mesh[pointer].coordZ<=SubstrateThickness+BOX ){
                //Burried
                if(mesh[pointer].coordZ>SubstrateThickness && mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz+Tox ){
                    material[pointer].Type=2;
                    material[pointer].k=SiO2_permi;
                    material[pointer].phi=volB;
                    material[pointer].dop=0;
                    material[pointer].u=0;
                    material[pointer].v=0;
                    material[pointer].mun=0;
                    material[pointer].mup=0;
                    material[pointer].tau=0;
                    material[pointer].r=0;
                    material[pointer].rho=0;
                }

                //NWR oxide
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz+Tox && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery1+NWRradiusy+Tox && mesh[pointer].coordY>=lx/2-NWRradiusy-Tox){
                        material[pointer].Type=2;
                        material[pointer].k=SiO2_permi;
                        material[pointer].phi=volB;
                        material[pointer].dop=0;
                        material[pointer].u=0;
                        material[pointer].v=0;
                        material[pointer].mun=0;
                        material[pointer].mup=0;
                        material[pointer].tau=0;
                        material[pointer].r=0;
                        material[pointer].rho=0;
                    }
                }

                //NWR1 channel
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery1+NWRradiusy && mesh[pointer].coordY>=NWRCentery1-NWRradiusy){
                        //P+
                        material[pointer].Type=1;
                        material[pointer].k=Si_permi;
                        material[pointer].dop=-Nai;
                        material[pointer].phi=((volD+volS)-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                        material[pointer].u=exp((-1)*(volD+volS)/VT);
                        material[pointer].v=exp((volD+volS)/VT);
                        material[pointer].mun=munCal(Tamb, Na, material[pointer].Type); // max Na Nd
                        material[pointer].mup=mupCal(Tamb, Na, material[pointer].Type);
                        material[pointer].tau=tauPCal(0);
                        material[pointer].r=SRHrecomb3D(i,j,k);
                        material[pointer].rho=0;

                        if(mesh[pointer].coordX<JunctionLength){
                            material[pointer].Type=1;
                            material[pointer].k=Si_permi;
                            material[pointer].dop=-NaPlusi;
                            material[pointer].phi=(volS-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            material[pointer].u=exp((-1)*volS/VT);
                            material[pointer].v=exp(volS/VT);
                            material[pointer].mun=munCal(Tamb, NaPlus, material[pointer].Type); // max Na Nd
                            material[pointer].mup=mupCal(Tamb, NaPlus, material[pointer].Type);
                            material[pointer].tau=tauPCal(0);
                            material[pointer].r=SRHrecomb3D(i,j,k);
                            material[pointer].rho=0;

                        }
                        if(mesh[pointer].coordX>ly-JunctionLength){
                            material[pointer].Type=1;
                            material[pointer].k=Si_permi;
                            material[pointer].dop=-NaPlusi;
                            material[pointer].phi=(volD-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            material[pointer].u=exp((-1)*volD/VT);
                            material[pointer].v=exp(volD/VT);
                            material[pointer].mun=munCal(Tamb, NaPlus, material[pointer].Type); // max Na Nd
                            material[pointer].mup=mupCal(Tamb, NaPlus, material[pointer].Type);
                            material[pointer].tau=tauPCal(0);
                            material[pointer].r=SRHrecomb3D(i,j,k);
                            material[pointer].rho=0;
                        }
                    }
                }
            }
        }
    }

}

void DDmodel::DDInitialGuess2NWR3D(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                //substrate
                if(mesh[pointer].coordZ<=SubstrateThickness ){
                    material[pointer].Type=4;
                    material[pointer].k=Si_permi;
                    material[pointer].dop=0;
                    material[pointer].phi=volB;
                    material[pointer].u=exp((-1)*volB/VT);
                    material[pointer].v=exp(volB/VT);
                    material[pointer].mun=munCal(Tamb, 0, material[pointer].Type); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, material[pointer].Type);
                    material[pointer].tau=tauPCal(0);
                    material[pointer].r=SRHrecomb3D(i,j,k);
                    material[pointer].rho=0;
                }

                //electrolyte
                if(mesh[pointer].coordZ>SubstrateThickness+BOX){
                    material[pointer].Type=3;
                    material[pointer].k=Water_permi;
                    //n+
                    //material[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2-VT*log(0.5*Ndi_SD+sqrt(pow(0.5*Ndi_SD,2)+1)))+volG;
                    //P+
                    //material[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2+VT*log(0.5*Nai_SD+sqrt(pow(0.5*Nai_SD,2)+1)))+volG;
                    //metal
                    //phi_new[pointer]=(Si_chi+(Si_Eg)/2)-(wfG)+volG;
                    //Intrinsic
                    material[pointer].phi=(Si_chi+(Si_Eg)/2)-(Si_chi+(Si_Eg)/2)+volG;
                    material[pointer].dop=0;
                    material[pointer].u=exp((-1)*volG/VT);
                    material[pointer].v=exp(volG/VT);
                    material[pointer].mun=munCal(Tamb, 0, material[pointer].Type); // max Na Nd
                    material[pointer].mup=mupCal(Tamb, 0, material[pointer].Type);
                    material[pointer].tau=tauNCal(Si_ni_cm);
                    material[pointer].r=0;
                    material[pointer].rho=0;
                }

                //box oxide
                //Burried
                if(mesh[pointer].coordZ>SubstrateThickness && mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz+Tox ){
                    material[pointer].Type=2;
                    material[pointer].k=SiO2_permi;
                    material[pointer].phi=volB;
                    material[pointer].dop=0;
                    material[pointer].u=0;
                    material[pointer].v=0;
                    material[pointer].mun=0;
                    material[pointer].mup=0;
                    material[pointer].tau=0;
                    material[pointer].r=0;
                    material[pointer].rho=0;
                }

                //NWR oxide
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz+Tox && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery1+NWRradiusy+Tox && mesh[pointer].coordY>=lx/2-NWRradiusy-Tox){
                        material[pointer].Type=2;
                        material[pointer].k=SiO2_permi;
                        material[pointer].phi=volB;
                        material[pointer].dop=0;
                        material[pointer].u=0;
                        material[pointer].v=0;
                        material[pointer].mun=0;
                        material[pointer].mup=0;
                        material[pointer].tau=0;
                        material[pointer].r=0;
                        material[pointer].rho=0;
                    }
                }

                //NWR1 channel
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery1+NWRradiusy && mesh[pointer].coordY>=NWRCentery1-NWRradiusy){
                        //P+
                        material[pointer].Type=1;
                        material[pointer].k=Si_permi;
                        material[pointer].dop=-Nai;
                        material[pointer].phi=((volD+volS)-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                        material[pointer].u=exp((-1)*(volD+volS)/VT);
                        material[pointer].v=exp((volD+volS)/VT);
                        material[pointer].mun=munCal(Tamb, Na, material[pointer].Type); // max Na Nd
                        material[pointer].mup=mupCal(Tamb, Na, material[pointer].Type);
                        material[pointer].tau=tauPCal(0);
                        material[pointer].r=SRHrecomb3D(i,j,k);
                        material[pointer].rho=0;

                        if(mesh[pointer].coordX<JunctionLength){
                            material[pointer].Type=1;
                            material[pointer].k=Si_permi;
                            material[pointer].dop=-NaPlusi;
                            material[pointer].phi=(volS-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            material[pointer].u=exp((-1)*volS/VT);
                            material[pointer].v=exp(volS/VT);
                            material[pointer].mun=munCal(Tamb, NaPlus, material[pointer].Type); // max Na Nd
                            material[pointer].mup=mupCal(Tamb, NaPlus, material[pointer].Type);
                            material[pointer].tau=tauPCal(0);
                            material[pointer].r=SRHrecomb3D(i,j,k);
                            material[pointer].rho=0;

                        }
                        if(mesh[pointer].coordX>ly-JunctionLength){
                            material[pointer].Type=1;
                            material[pointer].k=Si_permi;
                            material[pointer].dop=-NaPlusi;
                            material[pointer].phi=(volD-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            material[pointer].u=exp((-1)*volD/VT);
                            material[pointer].v=exp(volD/VT);
                            material[pointer].mun=munCal(Tamb, NaPlus, material[pointer].Type); // max Na Nd
                            material[pointer].mup=mupCal(Tamb, NaPlus, material[pointer].Type);
                            material[pointer].tau=tauPCal(0);
                            material[pointer].r=SRHrecomb3D(i,j,k);
                            material[pointer].rho=0;
                        }
                    }
                }

                //NWR2 channel
                if(mesh[pointer].coordZ<=SubstrateThickness+BOX+2*NWRradiusz && mesh[pointer].coordZ>=SubstrateThickness+BOX){
                    if(mesh[pointer].coordY<=NWRCentery2+NWRradiusy && mesh[pointer].coordY>=NWRCentery2-NWRradiusy){
                        //P+
                        material[pointer].Type=1;
                        material[pointer].k=Si_permi;
                        material[pointer].dop=-Nai;
                        material[pointer].phi=((volD+volS)-VT*log(0.5*Nai+sqrt(pow(0.5*Nai,2)+1)));
                        material[pointer].u=exp((-1)*(volD+volS)/VT);
                        material[pointer].v=exp((volD+volS)/VT);
                        material[pointer].mun=munCal(Tamb, Na, material[pointer].Type); // max Na Nd
                        material[pointer].mup=mupCal(Tamb, Na, material[pointer].Type);
                        material[pointer].tau=tauPCal(0);
                        material[pointer].r=SRHrecomb3D(i,j,k);
                        material[pointer].rho=0;

                        if(mesh[pointer].coordX<JunctionLength){
                            material[pointer].Type=1;
                            material[pointer].k=Si_permi;
                            material[pointer].dop=-NaPlusi;
                            material[pointer].phi=(volS-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            material[pointer].u=exp((-1)*volS/VT);
                            material[pointer].v=exp(volS/VT);
                            material[pointer].mun=munCal(Tamb, NaPlus, material[pointer].Type); // max Na Nd
                            material[pointer].mup=mupCal(Tamb, NaPlus, material[pointer].Type);
                            material[pointer].tau=tauPCal(0);
                            material[pointer].r=SRHrecomb3D(i,j,k);
                            material[pointer].rho=0;

                        }
                        if(mesh[pointer].coordX>ly-JunctionLength){
                            material[pointer].Type=1;
                            material[pointer].k=Si_permi;
                            material[pointer].dop=-NaPlusi;
                            material[pointer].phi=(volD-VT*log(0.5*NaPlusi+sqrt(pow(0.5*NaPlusi,2)+1)));
                            material[pointer].u=exp((-1)*volD/VT);
                            material[pointer].v=exp(volD/VT);
                            material[pointer].mun=munCal(Tamb, NaPlus, material[pointer].Type); // max Na Nd
                            material[pointer].mup=mupCal(Tamb, NaPlus, material[pointer].Type);
                            material[pointer].tau=tauPCal(0);
                            material[pointer].r=SRHrecomb3D(i,j,k);
                            material[pointer].rho=0;
                        }
                    }
                }
            }
        }
    }

}

double DDmodel::PoissonSolver2D(){

    DD_loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        DD_loop++;

        errPhi=PoissonGaussSeidel2D();

        if(errPhi_max < errPhi) {errPhi_max=errPhi;}

        if(DD_loop%1000==0)
        cout <<"PS:"<< DD_loop <<"\t" <<errPhi<<"\t"<<errPhi_max<<endl;

        if(DD_loop%100000==0)
        PrintMaterial2D("Poisson_temp.txt");

    }while(errPhi>SimTolPoisson);

    return errPhi_max;

}

double DDmodel::PoissonGaussSeidel2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double phik=material[pointer].phi;

            material[pointer].phi=PoissonGaussSeidelInner2D(i,j);

            material[pointer].r=SRHrecomb2D(i,j);

            double error=abs(material[pointer].phi-phik);

            error=error/(abs(phik)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    PoissonBC2D();

    return max_val;

}

double DDmodel::PoissonGaussSeidelInner2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double permitivity_ip=material[pointer].k*material[pointer_ip].k / (0.5*material[pointer_ip].k+0.5*material[pointer].k);
    double permitivity_in=material[pointer].k*material[pointer_in].k / (0.5*material[pointer_in].k+0.5*material[pointer].k);
    double permitivity_jp=material[pointer].k*material[pointer_jp].k / (0.5*material[pointer_jp].k+0.5*material[pointer].k);
    double permitivity_jn=material[pointer].k*material[pointer_jn].k / (0.5*material[pointer_jn].k+0.5*material[pointer].k);

    double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
    double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
    double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
    double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
    double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
    double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);

    double f,df,phik;
    double volume=deltax*deltay;

    phik=material[pointer].phi;

    int Temp=material[pointer].Type;

    switch(Temp){

    case 1: //1=channel
        f=volume*ni_nm*(material[pointer].u*exp(phik/VT)-material[pointer].v*exp((-1)*phik/VT)-material[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(material[pointer].u*exp(phik/VT)+material[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 2: //2=insulator
        f=volume*material[pointer].rho/e0;
        df=0;
        break;
    case 3: //3=electrolyte
        f=volume*C0*(material[pointer].u*exp(phik/VT)-material[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0;
       df=volume*C0*(material[pointer].u*exp(phik/VT)+material[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 4: //4=Substrate
        f=volume*ni_nm*(material[pointer].u*exp(phik/VT)-material[pointer].v*exp((-1)*phik/VT)-material[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(material[pointer].u*exp(phik/VT)+material[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 5: //5=Analyte
        f=volume*material[pointer].rho/e0;
        df=0;
        break;
    default:
        cout << "Error! Undefined Type for Poisson Solver @ PoissonGaussSeidelInner2D."<<endl;
        exit(0);
    }

    return (((permitivity_ip*material[pointer_ip].phi/xstep_p+permitivity_in*material[pointer_in].phi/xstep_n)*deltay
            +(permitivity_jp*material[pointer_jp].phi/ystep_p+permitivity_jn*material[pointer_jn].phi/ystep_n)*deltax
            + f - df*phik )
            /
            ((permitivity_ip/xstep_p+permitivity_in/xstep_n)*deltay
            +(permitivity_jp/ystep_p+permitivity_jn/ystep_n)*deltax - df));
}

double DDmodel::PoissonSolver3D(){

    DD_loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        DD_loop++;

        errPhi=PoissonGaussSeidel3D();

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

double DDmodel::PoissonGaussSeidel3D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<pz-1; k++) {

                int pointer = (py)*(px)*(k) + (px)*(j) + (i);

                double phik=material[pointer].phi;

                material[pointer].phi=PoissonGaussSeidelInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=abs(material[pointer].phi-phik);

                error=error/(abs(phik)+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    PoissonBC3D();

    return max_val;

}

double DDmodel::PoissonGaussSeidelInner3D(int i, int j, int k){

    int pointer = (py)*(px)*(k) + (px)*(j) + (i);
    int pointer_ip = (py)*(px)*(k) + (px)*(j) + (i+1);
    int pointer_in = (py)*(px)*(k) + (px)*(j) + (i-1);
    int pointer_jp = (py)*(px)*(k) + (px)*(j+1) + (i);
    int pointer_jn = (py)*(px)*(k) + (px)*(j-1) + (i);
    int pointer_kp = (py)*(px)*(k+1) + (px)*(j) + (i);
    int pointer_kn = (py)*(px)*(k-1) + (px)*(j) + (i);

    double permitivity_ip=material[pointer].k*material[pointer_ip].k / (0.5*material[pointer_ip].k+0.5*material[pointer].k);
    double permitivity_in=material[pointer].k*material[pointer_in].k / (0.5*material[pointer_in].k+0.5*material[pointer].k);
    double permitivity_jp=material[pointer].k*material[pointer_jp].k / (0.5*material[pointer_jp].k+0.5*material[pointer].k);
    double permitivity_jn=material[pointer].k*material[pointer_jn].k / (0.5*material[pointer_jn].k+0.5*material[pointer].k);
    double permitivity_kp=material[pointer].k*material[pointer_kp].k / (0.5*material[pointer_kp].k+0.5*material[pointer].k);
    double permitivity_kn=material[pointer].k*material[pointer_kn].k / (0.5*material[pointer_kn].k+0.5*material[pointer].k);

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

    phik=material[pointer].phi;

    int Temp=material[pointer].Type;

    switch(Temp){

    case 1: //1=channel
        f=volume*ni_nm*(material[pointer].u*exp(phik/VT)-material[pointer].v*exp((-1)*phik/VT)-material[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(material[pointer].u*exp(phik/VT)+material[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 2: //2=insulator
        f=volume*material[pointer].rho/e0;
        df=0;
        break;
    case 3: //3=electrolyte
        f=volume*C0*(material[pointer].u*exp(phik/VT)-material[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0;
       df=volume*C0*(material[pointer].u*exp(phik/VT)+material[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 4: //4=Substrate
        f=volume*ni_nm*(material[pointer].u*exp(phik/VT)-material[pointer].v*exp((-1)*phik/VT)-material[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(material[pointer].u*exp(phik/VT)+material[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
        break;
    case 5: //5=Analyte
        f=volume*material[pointer].rho/e0;
        df=0;
        break;
    default:
        cout << "Error! Undefined Type for Poisson Solver @ PoissonGaussSeidelInner3D."<<endl;
        exit(0);
    }

    return ((permitivity_ip*material[pointer_ip].phi/xstep_p+permitivity_in*material[pointer_in].phi/xstep_n)*deltay*deltaz
           +(permitivity_jp*material[pointer_jp].phi/ystep_p+permitivity_jn*material[pointer_jn].phi/ystep_n)*deltax*deltaz
           +(permitivity_kp*material[pointer_kp].phi/zstep_p+permitivity_kn*material[pointer_kn].phi/zstep_n)*deltax*deltay + f - df*phik )
            /
            ((permitivity_ip/xstep_p+permitivity_in/xstep_n)*deltay*deltaz
            +(permitivity_jp/ystep_p+permitivity_jn/ystep_n)*deltax*deltaz
            +(permitivity_kp/zstep_p+permitivity_kn/zstep_n)*deltax*deltay - df);
}

double DDmodel::ECSolver2D(){

    DD_loop=0;
    double errEC(0),errEC_max(0);

    do{
        DD_loop++;


        switch(StructureFlag){

        case 1:
            errEC=ECTypeA2D();
            break;
        case 2:
            errEC=ECTypeA2D();
            break;
        case 3:
            errEC=ECTypeB2D();
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

double DDmodel::ECTypeA2D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double uk = material[pointer].u;

            material[pointer].u=ECInner2D(i,j);

            material[pointer].r=SRHrecomb2D( i,j);

            double error=VT*abs(log(material[pointer].u)-log(uk));

            error=error/(abs(VT*log(uk))+1);

            if(error>max_val)
                max_val=error;
        }
    }

    ECBC2D();

    return max_val;

}

double DDmodel::ECTypeB2D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<SubstrateTOP; j++) {

            int pointer = (px)*(j) + (i);

            double uk = material[pointer].u;

            material[pointer].u=ECInner2D(i,j);

            material[pointer].r=SRHrecomb2D( i,j);

            double error=abs(log(material[pointer].u)-log(uk));

            if(error>max_val)
                max_val=error;
        }
    }

    /*
#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=ElectrolyteBottom+1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double uk = material[pointer].u;

            material[pointer].u=ECInner(i,j);

            material[pointer].r=SRHrecomb2D( i,j);

            double error=abs(log(material[pointer].u)-log(uk));

            if(error>max_val)
                max_val=error;
        }
    }
*/
    ECBC2D();

    return max_val;

}

double DDmodel::ECInner2D(int i, int j){

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

    double uf = material[pointer].u;

    double Bip = (material[pointer].mun+material[pointer_ip].mun)/2.0*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT);
    double Bin = (material[pointer].mun+material[pointer_in].mun)/2.0*Bern(material[pointer_in].phi/VT-material[pointer].phi/VT, material[pointer_in].phi/VT);
    double Bjp = (material[pointer].mun+material[pointer_jp].mun)/2.0*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, material[pointer_jp].phi/VT);
    double Bjn = (material[pointer].mun+material[pointer_jn].mun)/2.0*Bern(material[pointer_jn].phi/VT-material[pointer].phi/VT, material[pointer_jn].phi/VT);

    double f=volume/VT*material[pointer].r;

    double df=volume/VT*(pow(material[pointer].v,2)*exp(-material[pointer].phi/VT)+2*material[pointer].v+exp(material[pointer].phi/VT))/material[pointer].tau/
                pow(material[pointer].u*exp(material[pointer].phi/VT)+material[pointer].v*exp(-material[pointer].phi/VT)+2,2);

    uf=((Bip*material[pointer_ip].u/xstep_p+Bin*material[pointer_in].u/xstep_n)*deltay
       +(Bjp*material[pointer_jp].u/ystep_p+Bjn*material[pointer_jn].u/ystep_n)*deltax
       -f+df*uf)/((Bip/xstep_p+Bin/xstep_n)*deltay+(Bjp/ystep_p+Bjn/ystep_n)*deltax+df);

    return uf;
}

double DDmodel::HCSolver2D(){

    DD_loop=0;
    double errHC(0),errHC_max(0);

    do{
        DD_loop++;

        switch(StructureFlag){

        case 1:
            errHC=HCTypeA2D();
            break;
        case 2:
            errHC=HCTypeA2D();
            break;
        case 3:
            errHC=HCTypeB2D();
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

double DDmodel::HCTypeA2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double vk = material[pointer].v;

            material[pointer].v=HCInner2D(i,j);

            material[pointer].r=SRHrecomb2D( i,j);

            double error=VT*abs(log(material[pointer].v)-log(vk));

            error=error/(abs(VT*log(vk))+1);

            if(error>max_val)
                max_val=error;
        }
    }

    HCBC2D();

    return max_val;
}

double DDmodel::HCTypeB2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<SubstrateTOP; j++) {

            int pointer = (px)*(j) + (i);

            double vk = material[pointer].v;

            material[pointer].v=HCInner2D(i,j);

            material[pointer].r=SRHrecomb2D( i,j);

            double error=abs(log(material[pointer].v)-log(vk));

            if(error>max_val)
                max_val=error;
        }
    }

    /*
#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=ElectrolyteBottom+1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double vk = material[pointer].v;

            material[pointer].v=HCInner(i,j);

            material[pointer].r=SRHrecomb2D( i,j);

            double error=abs(log(material[pointer].v)-log(vk));

            if(error>max_val)
                max_val=error;
        }
    }
    */
    HCBC2D();

    return max_val;
}

double DDmodel::HCInner2D(int i, int j){

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

    double vf = material[pointer].v;

    double Bip = (material[pointer].mup+material[pointer_ip].mup)/2.0*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);
    double Bin = (material[pointer].mup+material[pointer_in].mup)/2.0*Bern(material[pointer_in].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);
    double Bjp = (material[pointer].mup+material[pointer_jp].mup)/2.0*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);
    double Bjn = (material[pointer].mup+material[pointer_jn].mup)/2.0*Bern(material[pointer_jn].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);

    double f=volume/VT*material[pointer].r;

    double df=volume/VT*(pow(material[pointer].u,2)*exp(material[pointer].phi/VT)+2*material[pointer].u+exp(-material[pointer].phi/VT))/material[pointer].tau/
                pow(material[pointer].u*exp(material[pointer].phi/VT)+material[pointer].v*exp(-material[pointer].phi/VT)+2,2);

    vf =((Bip*material[pointer_ip].v/xstep_p+Bin*material[pointer_in].v/xstep_n)*deltay
         +(Bjp*material[pointer_jp].v/ystep_p+Bjn*material[pointer_jn].v/ystep_n)*deltax
         -f+df*vf)/((Bip/xstep_p+Bin/xstep_n)*deltay+(Bjp/ystep_p+Bjn/ystep_n)*deltax+df);

    return vf;
}

void DDmodel::PoissonBC2D(){

    switch(StructureFlag){

    case 1:
        PoissonBC2D_PN();
        break;
    case 2:
        PoissonBC2D_MOSFET();
        break;
    case 3:
        PoissonBC2D_ISFET();
        break;
    default:
        cout <<"Undefined Boundary @ PoissonBC2D."<<endl;
        exit(0);
    }
}

void DDmodel::PoissonBC2D_PN(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);
            material[pointer1].phi=material[pointer2].phi;

            pointer1 = (px)*(py-1) + (i);
            pointer2 = (px)*(py-2) + (i);
            material[pointer1].phi=material[pointer2].phi;
        }
}

void DDmodel::PoissonBC2D_MOSFET(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);

            if( mesh[pointer1].coordX > JunctionLength && mesh[pointer1].coordX < lx-JunctionLength ){

                double ystep=abs(mesh[pointer1].coordY-mesh[pointer2].coordY);
                double qfactor=Si_permi/SiO2_permi*Tox/ystep;

                material[pointer1].phi=(volG+qfactor*material[pointer2].phi)/(1.0+qfactor);
            }
        }

#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);
            material[pointer1].phi=material[pointer2].phi;

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);
            material[pointer1].phi=material[pointer2].phi;
        }
}

void DDmodel::PoissonBC2D_ISFET(){

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

                material[pointer1].phi=material[pointer2].phi;

            }

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                material[pointer1].phi=material[pointer2].phi;
            }
        }
}

void DDmodel::ECBC2D(){

    switch(StructureFlag){

    case 1:
        ECBC2D_PN();
        break;
    case 2:
        ECBC2D_MOSFET();
        break;
    case 3:
        ECBC2D_ISFET();
        break;
    default:
        cout <<"Undefined Boundary @ ECBC2D."<<endl;
        exit(0);
    }
}

void DDmodel::HCBC2D(){

    switch(StructureFlag){

    case 1:
        HCBC2D_PN();
        break;
    case 2:
        HCBC2D_MOSFET();
        break;
    case 3:
        HCBC2D_ISFET();
        break;
    default:
        cout <<"Undefined Boundary @ HCBC2D."<<endl;
        exit(0);
    }
}

void DDmodel::ECBC2D_PN(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);
            material[pointer1].u=material[pointer2].u;
            material[pointer1].r=material[pointer2].r;

            pointer1 = (px)*(py-1) + (i);
            pointer2 = (px)*(py-2) + (i);
            material[pointer1].u=material[pointer2].u;
            material[pointer1].r=material[pointer2].r;
        }
}

void DDmodel::HCBC2D_PN(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);
            material[pointer1].v=material[pointer2].v;
            material[pointer1].r=material[pointer2].r;

            pointer1 = (px)*(py-1) + (i);
            pointer2 = (px)*(py-2) + (i);
            material[pointer1].v=material[pointer2].v;
            material[pointer1].r=material[pointer2].r;
        }
}

void DDmodel::ECBC2D_MOSFET(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);

            if( mesh[pointer1].coordX > JunctionLength && mesh[pointer1].coordX < lx-JunctionLength ){

                material[pointer1].u=material[pointer2].u;
                material[pointer1].r=material[pointer2].r;
            }
        }

#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);
            material[pointer1].u=material[pointer2].u;
            material[pointer1].r=material[pointer2].r;

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);
            material[pointer1].u=material[pointer2].u;
            material[pointer1].r=material[pointer2].r;
        }
}

void DDmodel::HCBC2D_MOSFET(){

#pragma omp parallel for
        for (int i=0; i<px; i++) {

            int pointer1 = (px)*(0) + (i);
            int pointer2 = (px)*(1) + (i);

            if( mesh[pointer1].coordX > JunctionLength && mesh[pointer1].coordX < lx-JunctionLength ){

                material[pointer1].v=material[pointer2].v;
                material[pointer1].r=material[pointer2].r;
            }
        }

#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);
            material[pointer1].v=material[pointer2].v;
            material[pointer1].r=material[pointer2].r;

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);
            material[pointer1].v=material[pointer2].v;
            material[pointer1].r=material[pointer2].r;
        }
}

void DDmodel::ECBC2D_ISFET(){

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
            material[pointer1].u=material[pointer2].u;
            material[pointer1].r=material[pointer2].r;

            pointer1 = (px)*(ElectrolyteBottom) + (i);
            pointer2 = (px)*(ElectrolyteBottom+1) + (i);
            material[pointer1].u=material[pointer2].u;
            material[pointer1].r=material[pointer2].r;
        }
#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                material[pointer1].u=material[pointer2].u;
                material[pointer1].r=material[pointer2].r;

            }

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                material[pointer1].u=material[pointer2].u;
                material[pointer1].r=material[pointer2].r;
            }
        }
}

void DDmodel::HCBC2D_ISFET(){

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
            material[pointer1].v=material[pointer2].v;
            material[pointer1].r=material[pointer2].r;

            pointer1 = (px)*(ElectrolyteBottom) + (i);
            pointer2 = (px)*(ElectrolyteBottom+1) + (i);
            material[pointer1].v=material[pointer2].v;
            material[pointer1].r=material[pointer2].r;
        }
#pragma omp parallel for
        for (int j=0; j<py; j++) {

            int pointer1 = (px)*(j) + (0);
            int pointer2 = (px)*(j) + (1);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                material[pointer1].v=material[pointer2].v;
                material[pointer1].r=material[pointer2].r;

            }

            pointer1 = (px)*(j) + (px-1);
            pointer2 = (px)*(j) + (px-2);

            if(mesh[pointer1].coordY<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordY>SubstrateThickness){

                material[pointer1].v=material[pointer2].v;
                material[pointer1].r=material[pointer2].r;
            }
        }
}

void DDmodel::ECBC3D(){

    switch(StructureFlag){

    case 1:
        ECBC3D_PN();
        break;
    case 2:
        ECBC3D_MOSFET();
        break;
    case 3:
        ECBC3D_ISFET();
        break;
    case 4:
        ECBC3D_1NWR();
        break;
    case 5:
        ECBC3D_2NWR();
        break;
    default:
        cout <<"Undefined Boundary @ ECBC3D."<<endl;
        exit(0);
    }
}

void DDmodel::ECBC3D_PN(){

    for(int i=0;i<px;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(0) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(1) +  (px)*(j) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(pz-1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(pz-2) +  (px)*(j) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(0) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(py-1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(py-2) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }
}

void DDmodel::ECBC3D_MOSFET(){

    for(int i=0;i<px;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(0) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(1) +  (px)*(j) + (i);

            if(mesh[pointer0].coordX > JunctionLength && mesh[pointer0].coordX < lx-JunctionLength ){
                material[pointer0].u=material[pointer1].u;
                material[pointer0].r=material[pointer1].r;
            }
        }
    }

    for(int j=0;j<py;j++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(j) + (0);
            int pointer1 =(px)*(py)*(k) +  (px)*(j) + (1);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(pz-1) +  (px)*(j) + (px-1);
            pointer1 =(px)*(py)*(pz-2) +  (px)*(j) + (px-2);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(0) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(py-1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(py-2) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }
}

void DDmodel::ECBC3D_ISFET(){

    for (int i=0; i<px; i++) {
        for (int j=0;j<py;j++){

            int pointer1 = (py)*(px)*(SubstrateTOP) + (px)*(j) + (i);
            int pointer2 = (py)*(px)*(SubstrateTOP-1) + (px)*(j) + (i);
            material[pointer1].u=material[pointer2].u;
            material[pointer1].r=material[pointer2].r;

            pointer1 = (py)*(px)*(ElectrolyteBottom) + (px)*(j) + (i);
            pointer2 = (py)*(px)*(ElectrolyteBottom+1) + (px)*(j) + (i);
            material[pointer1].u=material[pointer2].u;
            material[pointer1].r=material[pointer2].r;
        }
    }

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(mesh[pointer0].coordZ<=SubstrateThickness-JunctionDepth || mesh[pointer0].coordZ>SubstrateThickness ){
                material[pointer0].u=material[pointer1].u;
                material[pointer0].r=material[pointer1].r;
            }


            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(mesh[pointer0].coordZ<=SubstrateThickness-JunctionDepth || mesh[pointer0].coordZ>SubstrateThickness ){
                material[pointer0].u=material[pointer1].u;
                material[pointer0].r=material[pointer1].r;
            }
        }
    }
}

void DDmodel::ECBC3D_1NWR(){

    for(int i=0;i<px;i++){
        for(int j=NWRleft1;j<=NWRright1;j++){
            int pointer0 =(px)*(py)*(NWRbottom1) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom1+1) +  (px)*(j) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop1-1) +  (px)*(j) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=NWRbottom1;k<=NWRtop1;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(NWRleft1) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(NWRleft1+1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(NWRright1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(NWRright1-1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }
}

void DDmodel::ECBC3D_2NWR(){

    //NWR1
    for(int i=0;i<px;i++){
        for(int j=NWRleft1;j<=NWRright1;j++){
            int pointer0 =(px)*(py)*(NWRbottom1) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom1+1) +  (px)*(j) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop1-1) +  (px)*(j) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=NWRbottom1;k<=NWRtop1;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(NWRleft1) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(NWRleft1+1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(NWRright1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(NWRright1-1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }

    //NWR2
    for(int i=0;i<px;i++){
        for(int j=NWRleft2;j<=NWRright2;j++){
            int pointer0 =(px)*(py)*(NWRbottom2) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom2+1) +  (px)*(j) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop2) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop2-1) +  (px)*(j) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=NWRbottom2;k<=NWRtop2;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(NWRleft2) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(NWRleft2+1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(NWRright2) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(NWRright2-1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;
        }
    }
}

void DDmodel::HCBC3D(){

    switch(StructureFlag){

    case 1:
        HCBC3D_PN();
        break;
    case 2:
        HCBC3D_MOSFET();
        break;
    case 3:
        HCBC3D_ISFET();
        break;
    case 4:
        HCBC3D_1NWR();
        break;
    case 5:
        HCBC3D_2NWR();
        break;
    default:
        cout <<"Undefined Boundary @ ECBC2D."<<endl;
        exit(0);
    }
}

void DDmodel::HCBC3D_PN(){

    for(int i=0;i<px;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(0) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(1) +  (px)*(j) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(pz-1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(pz-2) +  (px)*(j) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(0) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(1) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(py-1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(py-2) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }
}

void DDmodel::HCBC3D_MOSFET(){

    for(int i=0;i<px;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(0) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(1) +  (px)*(j) + (i);

            if(mesh[pointer0].coordX > JunctionLength && mesh[pointer0].coordX < lx-JunctionLength ){
                material[pointer0].v=material[pointer1].v;
                material[pointer0].r=material[pointer1].r;
            }
        }
    }

    for(int j=0;j<py;j++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(j) + (0);
            int pointer1 =(px)*(py)*(k) +  (px)*(j) + (1);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(pz-1) +  (px)*(j) + (px-1);
            pointer1 =(px)*(py)*(pz-2) +  (px)*(j) + (px-2);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=0;k<pz;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(0) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(1) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(py-1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(py-2) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }
}

void DDmodel::HCBC3D_ISFET(){
    for (int i=0; i<px; i++) {
        for (int j=0;j<py;j++){

            int pointer1 = (py)*(px)*(SubstrateTOP) + (px)*(j) + (i);
            int pointer2 = (py)*(px)*(SubstrateTOP-1) + (px)*(j) + (i);
            material[pointer1].v=material[pointer2].v;
            material[pointer1].r=material[pointer2].r;

            pointer1 = (py)*(px)*(ElectrolyteBottom) + (px)*(j) + (i);
            pointer2 = (py)*(px)*(ElectrolyteBottom+1) + (px)*(j) + (i);
            material[pointer1].v=material[pointer2].v;
            material[pointer1].r=material[pointer2].r;
        }
    }

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            material[pointer0].u=material[pointer1].u;
            material[pointer0].r=material[pointer1].r;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(mesh[pointer1].coordZ<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordZ>SubstrateThickness){

                material[pointer0].v=material[pointer1].v;
                material[pointer0].r=material[pointer1].r;
            }


            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(mesh[pointer1].coordZ<=SubstrateThickness-JunctionDepth|| mesh[pointer1].coordZ>SubstrateThickness){

                material[pointer0].v=material[pointer1].v;
                material[pointer0].r=material[pointer1].r;
            }
        }
    }
}

void DDmodel::HCBC3D_1NWR(){

    for(int i=0;i<px;i++){
        for(int j=NWRleft1;j<=NWRright1;j++){
            int pointer0 =(px)*(py)*(NWRbottom1) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom1+1) +  (px)*(j) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop1-1) +  (px)*(j) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int i=0;i<px;i++){
        for(int k=NWRbottom1;k<=NWRtop1;k++){
            int pointer0 =(px)*(py)*(k) +  (px)*(NWRleft1) + (i);
            int pointer1 =(px)*(py)*(k) +  (px)*(NWRleft1+1) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(NWRright1) + (i);
            pointer1 =(px)*(py)*(k) +  (px)*(NWRright1-1) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }
}

void DDmodel::HCBC3D_2NWR(){

    //NWR1
    for(int i=NWRleft1;i<=NWRright1;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(NWRbottom1) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom1+1) +  (px)*(j) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop1) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop1-1) +  (px)*(j) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int k=NWRbottom1;k<=NWRtop1;k++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(k) +  (px)*(j) + (NWRleft1);
            int pointer1 =(px)*(py)*(k) +  (px)*(j) + (NWRleft1+1);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(j) + (NWRright1);
            pointer1 =(px)*(py)*(k) +  (px)*(j) + (NWRright1-1);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }

    //NWR2
    for(int i=NWRleft2;i<=NWRright2;i++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(NWRbottom2) +  (px)*(j) + (i);
            int pointer1 =(px)*(py)*(NWRbottom2+1) +  (px)*(j) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(NWRtop2) +  (px)*(j) + (i);
            pointer1 =(px)*(py)*(NWRtop2-1) +  (px)*(j) + (i);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }

    for(int k=NWRbottom2;k<=NWRtop2;k++){
        for(int j=0;j<py;j++){
            int pointer0 =(px)*(py)*(k) +  (px)*(j) + (NWRleft2);
            int pointer1 =(px)*(py)*(k) +  (px)*(j) + (NWRleft2+1);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;

            pointer0 =(px)*(py)*(k) +  (px)*(j) + (NWRright2);
            pointer1 =(px)*(py)*(k) +  (px)*(j) + (NWRright2-1);

            material[pointer0].v=material[pointer1].v;
            material[pointer0].r=material[pointer1].r;
        }
    }
}

void DDmodel::PoissonBC3D(){

    switch(StructureFlag){

    case 1:
        PoissonBC3D_PN();
        break;
    case 2:
        PoissonBC3D_MOSFET();
        break;
    case 3:
        PoissonBC3D_ISFET();
        break;
    case 4:
        PoissonBC3D_1NWR();
        break;
    case 5:
        PoissonBC3D_2NWR();
        break;
    default:
        cout <<"Undefined Boundary @ PoissonBC3D."<<endl;
        exit(0);
    }

}

void DDmodel::PoissonBC3D_PN(){

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer0 = (py)*(px)*(0) + (px)*(j) + (i);
            int pointer1 = (py)*(px)*(1) + (px)*(j) + (i);
            material[pointer0].phi=material[pointer1].phi;

            pointer0 = (py)*(px)*(pz-1) + (px)*(j) + (i);
            pointer1 = (py)*(px)*(pz-2) + (px)*(j) + (i);
            material[pointer0].phi=material[pointer1].phi;
        }
    }

#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);
            material[pointer0].phi=material[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);
            material[pointer0].phi=material[pointer1].phi;
        }
    }
}

void DDmodel::PoissonBC3D_MOSFET(){


#pragma omp parallel for
        for (int i=0; i<px; i++) {
            for (int j=0;j<py;j++){

                int pointer1 = (py)*(px)*(0) + (px)*(j) + (i);
                int pointer2 = (py)*(px)*(1) + (px)*(j) + (i);

                if(mesh[pointer1].coordX > JunctionLength && mesh[pointer1].coordX < lx-JunctionLength ){

                    double zstep=abs(mesh[pointer1].coordZ-mesh[pointer2].coordZ);
                    double qfactor=Si_permi/SiO2_permi*Tox/zstep;

                    material[pointer1].phi=(volG+qfactor*material[pointer2].phi)/(1.0+qfactor);
                }
            }
        }

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);
            material[pointer0].phi=material[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);
            material[pointer0].phi=material[pointer1].phi;
        }
    }

    for (int k=0;k<pz;k++){
        for (int j=0;j<py;j++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);
            material[pointer0].phi=material[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);
            material[pointer0].phi=material[pointer1].phi;
        }
    }
}

void DDmodel::PoissonBC3D_ISFET(){

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            material[pointer0].phi=material[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            material[pointer0].phi=material[pointer1].phi;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(mesh[pointer0].coordZ<=SubstrateThickness-JunctionDepth || mesh[pointer0].coordZ>SubstrateThickness ){
                material[pointer0].phi=material[pointer1].phi;
            }

            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(mesh[pointer0].coordZ<=SubstrateThickness-JunctionDepth || mesh[pointer0].coordZ>SubstrateThickness ){
                material[pointer0].phi=material[pointer1].phi;
            }
        }
    }
}

void DDmodel::PoissonBC3D_1NWR(){

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            material[pointer0].phi=material[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            material[pointer0].phi=material[pointer1].phi;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(material[pointer0].Type!=1){
                material[pointer0].phi=material[pointer1].phi;
            }

            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(material[pointer0].Type!=1){
                material[pointer0].phi=material[pointer1].phi;
            }
        }
    }
}

void DDmodel::PoissonBC3D_2NWR(){

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer1 = (py)*(px)*(k) + (px)*(1) + (i);

            material[pointer0].phi=material[pointer1].phi;

            pointer0 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer1 = (py)*(px)*(k) + (px)*(py-2) + (i);

            material[pointer0].phi=material[pointer1].phi;
        }
    }

    for (int j=0;j<py;j++){
        for (int k=0;k<pz;k++){

            int pointer0 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer1 = (py)*(px)*(k) + (px)*(j) + (1);

            if(material[pointer0].Type!=1){
                material[pointer0].phi=material[pointer1].phi;
            }

            pointer0 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-2);

            if(material[pointer0].Type!=1){
                material[pointer0].phi=material[pointer1].phi;
            }
        }
    }
}

void DDmodel::BernoulliX(){

    double x(1.0),v1,v2;

    do{
        x=x/2.0;
        v1=x/(exp(x)-1.0);
        v2=1.0-x/2.0;

    }while((v1-v2)>1e-10);
    bernXl=x;
}

void DDmodel::PrintMaterial2D(string path){

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
                   << material[pointer].k << '\t' <<material[pointer].dop << '\t' <<material[pointer].phi << '\t'
                   << material[pointer].u << '\t' << material[pointer].v << '\t' << material[pointer].r << '\t'
                   << material[pointer].rho << '\t'<< material[pointer].mun << '\t'<< material[pointer].mup << '\t'
                   << material[pointer].tau << '\t'<< material[pointer].Ex << '\t'<< material[pointer].Ey << '\t'
                   << material[pointer].Type << '\t'<< material[pointer].Crho <<endl;

        }
    }

    output.close();
}

void DDmodel::PrintMaterial3D(string path){

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
                       << material[pointer].k << '\t' <<material[pointer].dop << '\t' <<material[pointer].phi << '\t'
                       << material[pointer].u << '\t' << material[pointer].v << '\t' << material[pointer].r << '\t'
                       << material[pointer].rho << '\t'<< material[pointer].mun << '\t'<< material[pointer].mup << '\t'
                       << material[pointer].tau << '\t'<< material[pointer].Ex << '\t'<< material[pointer].Ey << '\t'
                       << material[pointer].Ez << '\t'<< material[pointer].Type << '\t'<< material[pointer].Crho <<endl;
            }
        }
    }

    output.close();
}

void DDmodel::ReadMaterial2D(string path){

    DDInitialize();

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
                  >> material[pointer].k >> material[pointer].dop >> material[pointer].phi
                  >> material[pointer].u >> material[pointer].v >> material[pointer].r
                  >> material[pointer].rho >> material[pointer].mun >> material[pointer].mup
                  >> material[pointer].tau >> material[pointer].Ex >> material[pointer].Ey
                  >> material[pointer].Type >> material[pointer].Crho;
        }
    }

    input.close();
}

void DDmodel::ReadMaterial3D(string path){

    DDInitialize();

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
                  >> material[pointer].k >> material[pointer].dop >> material[pointer].phi
                  >> material[pointer].u >> material[pointer].v >> material[pointer].r
                  >> material[pointer].rho >> material[pointer].mun >> material[pointer].mup
                  >> material[pointer].tau >> material[pointer].Ex >> material[pointer].Ey
                  >> material[pointer].Ez >> material[pointer].Type >> material[pointer].Crho;
        }
    }

    input.close();
}

void DDmodel::RhoCalculation2D(){

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

            double permitivity_ip=material[pointer].k*material[pointer_ip].k / (0.5*material[pointer_ip].k+0.5*material[pointer].k);
            double permitivity_in=material[pointer].k*material[pointer_in].k / (0.5*material[pointer_in].k+0.5*material[pointer].k);
            double permitivity_jp=material[pointer].k*material[pointer_jp].k / (0.5*material[pointer_jp].k+0.5*material[pointer].k);
            double permitivity_jn=material[pointer].k*material[pointer_jn].k / (0.5*material[pointer_jn].k+0.5*material[pointer].k);

            double charge_totalemp = (-1)*((permitivity_ip*material[pointer_ip].phi - permitivity_ip*material[pointer].phi)/xstep_p*deltay
                                          +(permitivity_in*material[pointer_in].phi - permitivity_in*material[pointer].phi)/xstep_n*deltay
                                          +(permitivity_jp*material[pointer_jp].phi - permitivity_jp*material[pointer].phi)/ystep_p*deltax
                                          +(permitivity_jn*material[pointer_jn].phi - permitivity_jn*material[pointer].phi)/ystep_n*deltax)
                                          /(deltax*deltay);

            material[pointer].Crho=charge_totalemp*e0;
        }
    }
}

void DDmodel::RhoCalculation3D(){

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

                double permitivity_ip=material[pointer].k*material[pointer_ip].k / (0.5*material[pointer_ip].k+0.5*material[pointer].k);
                double permitivity_in=material[pointer].k*material[pointer_in].k / (0.5*material[pointer_in].k+0.5*material[pointer].k);
                double permitivity_jp=material[pointer].k*material[pointer_jp].k / (0.5*material[pointer_jp].k+0.5*material[pointer].k);
                double permitivity_jn=material[pointer].k*material[pointer_jn].k / (0.5*material[pointer_jn].k+0.5*material[pointer].k);
                double permitivity_kp=material[pointer].k*material[pointer_kp].k / (0.5*material[pointer_kp].k+0.5*material[pointer].k);
                double permitivity_kn=material[pointer].k*material[pointer_kn].k / (0.5*material[pointer_kn].k+0.5*material[pointer].k);

                double deltax=abs(mesh[pointer_ip].coordX-mesh[pointer_in].coordX)/2;
                double deltay=abs(mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)/2;
                double deltaz=abs(mesh[pointer_kp].coordZ-mesh[pointer_kn].coordZ)/2;
                double xstep_p=abs(mesh[pointer_ip].coordX-mesh[pointer].coordX);
                double xstep_n=abs(mesh[pointer_in].coordX-mesh[pointer].coordX);
                double ystep_p=abs(mesh[pointer_jp].coordY-mesh[pointer].coordY);
                double ystep_n=abs(mesh[pointer_jn].coordY-mesh[pointer].coordY);
                double zstep_p=abs(mesh[pointer_kp].coordZ-mesh[pointer].coordZ);
                double zstep_n=abs(mesh[pointer_kn].coordZ-mesh[pointer].coordZ);

                double charge_totalemp = (-1)*((permitivity_ip*material[pointer_ip].phi - permitivity_ip*material[pointer].phi)/xstep_p*deltay*deltaz
                                              +(permitivity_in*material[pointer_in].phi - permitivity_in*material[pointer].phi)/xstep_n*deltay*deltaz
                                              +(permitivity_jp*material[pointer_jp].phi - permitivity_jp*material[pointer].phi)/ystep_p*deltax*deltaz
                                              +(permitivity_jn*material[pointer_jn].phi - permitivity_jn*material[pointer].phi)/ystep_n*deltax*deltaz
                                              +(permitivity_kp*material[pointer_kp].phi - permitivity_kp*material[pointer].phi)/zstep_p*deltax*deltay
                                              +(permitivity_kn*material[pointer_kn].phi - permitivity_kn*material[pointer].phi)/zstep_n*deltax*deltay)
                                              /(deltax*deltay*deltaz);

                material[pointer].Crho=charge_totalemp*e0;
                //Crho = [C/nm3]
            }
        }
    }
}

void DDmodel::EfieldCalculation2D(){

    #pragma omp parallel for
    for (int i=1;i<px-1;i++){
        for (int j=1;j<py-1;j++){

            int pointer = (px)*(j) + (i);
            int pointer_ip =   (px)*(j) + (i+1);
            int pointer_in =   (px)*(j) + (i-1);
            int pointer_jp =   (px)*(j+1) + (i);
            int pointer_jn =   (px)*(j-1) + (i);

            double Efx=(-1)*(material[pointer_ip].phi-material[pointer_in].phi)/((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)*1e-9);
            double Efy=(-1)*(material[pointer_jp].phi-material[pointer_jn].phi)/((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)*1e-9);

            material[pointer].Ex=Efx;
            material[pointer].Ey=Efy;
        }
    }

    for (int i=0;i<px;i++){

        int pointer1 = (px)*(0) + (i);
        int pointer2 = (px)*(1) + (i);
        material[pointer1].Ex=material[pointer2].Ex;
        material[pointer1].Ey=material[pointer2].Ey;

        pointer1 = (px)*(py-1) + (i);
        pointer2 = (px)*(py-2) + (i);
        material[pointer1].Ex=material[pointer2].Ex;
        material[pointer1].Ey=material[pointer2].Ey;
    }

    for (int j=0;j<py;j++){

        int pointer1 = (px)*(j) + (0);
        int pointer2 = (px)*(j) + (1);
        material[pointer1].Ex=material[pointer2].Ex;
        material[pointer1].Ey=material[pointer2].Ey;

        pointer1 = (px)*(j) + (px-1);
        pointer2 = (px)*(j) + (px-2);
        material[pointer1].Ex=material[pointer2].Ex;
        material[pointer1].Ey=material[pointer2].Ey;
    }

}

void DDmodel::EfieldCalculation3D(){

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

                double Efx=(-1)*(material[pointer_ip].phi-material[pointer_in].phi)/((mesh[pointer_ip].coordX-mesh[pointer_in].coordX)*1e-9);
                double Efy=(-1)*(material[pointer_jp].phi-material[pointer_jn].phi)/((mesh[pointer_jp].coordY-mesh[pointer_jn].coordY)*1e-9);
                double Efz=(-1)*(material[pointer_kp].phi-material[pointer_kn].phi)/((mesh[pointer_kp].coordZ-mesh[pointer_kn].coordZ)*1e-9);

                material[pointer].Ex=Efx;
                material[pointer].Ey=Efy;
                material[pointer].Ez=Efz;
            }
        }
    }

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer1 = (py)*(px)*(0) + (px)*(j) + (i);
            int pointer2 = (py)*(px)*(1) + (px)*(j) + (i);
            material[pointer1].Ex=material[pointer2].Ex;
            material[pointer1].Ey=material[pointer2].Ey;
            material[pointer1].Ez=material[pointer2].Ez;

            pointer1 = (py)*(px)*(pz-1) + (px)*(j) + (i);
            pointer2 = (py)*(px)*(pz-2) + (px)*(j) + (i);
            material[pointer1].Ex=material[pointer2].Ex;
            material[pointer1].Ey=material[pointer2].Ey;
            material[pointer1].Ez=material[pointer2].Ez;
        }
    }

    for (int i=0;i<px;i++){
        for (int k=0;k<pz;k++){

            int pointer1 = (py)*(px)*(k) + (px)*(0) + (i);
            int pointer2 = (py)*(px)*(k) + (px)*(1) + (i);
            material[pointer1].Ex=material[pointer2].Ex;
            material[pointer1].Ey=material[pointer2].Ey;
            material[pointer1].Ez=material[pointer2].Ez;

            pointer1 = (py)*(px)*(k) + (px)*(py-1) + (i);
            pointer2 = (py)*(px)*(k) + (px)*(py-2) + (i);
            material[pointer1].Ex=material[pointer2].Ex;
            material[pointer1].Ey=material[pointer2].Ey;
            material[pointer1].Ez=material[pointer2].Ez;
        }
    }

    for (int k=0;k<pz;k++){
        for (int j=0;j<py;j++){

            int pointer1 = (py)*(px)*(k) + (px)*(j) + (0);
            int pointer2 = (py)*(px)*(k) + (px)*(j) + (1);
            material[pointer1].Ex=material[pointer2].Ex;
            material[pointer1].Ey=material[pointer2].Ey;
            material[pointer1].Ez=material[pointer2].Ez;

            pointer1 = (py)*(px)*(k) + (px)*(j) + (px-1);
            pointer2 = (py)*(px)*(k) + (px)*(j) + (px-2);
            material[pointer1].Ex=material[pointer2].Ex;
            material[pointer1].Ey=material[pointer2].Ey;
            material[pointer1].Ez=material[pointer2].Ez;
        }
    }
}

void DDmodel::DDInitialize(){

    #pragma omp parallel for
    for(int i=0;i<L;i++){
        material[i].Crho=0;
        material[i].dop=0;
        material[i].Ex=0;
        material[i].Ey=0;
        material[i].Ez=0;
        material[i].k=0;
        material[i].mun=0;
        material[i].mup=0;
        material[i].phi=0;
        material[i].r=0;
        material[i].rho=0;
        material[i].tau=0;
        material[i].Type=0;
        material[i].u=0;
        material[i].v=0;
    }
}


double DDmodel::munCal(double T, double dopping, int f){

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

double DDmodel::mupCal(double T,double dopping, int f){

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

double DDmodel::tauNCal(double dopping){

    //Nd dopping
    //calculate hole lifetime
    return 3.94e-4/(1+dopping/1e21);
}

double DDmodel::tauPCal(double dopping){

    //Na dopping
    //calculate electron lifetime
    return 3.94e-4/(1+dopping/1e21);
}

double DDmodel::SRHrecomb2D(int i, int j){

    int pointer = (px)*(j) + (i);
    //http://www.iue.tuwien.ac.at/phd/entner/node11.html
    return (material[pointer].u*material[pointer].v-1)/material[pointer].tau/(material[pointer].u*exp(material[pointer].phi/VT)+material[pointer].v*exp((-1)*material[pointer].phi/VT)+2);
}

double DDmodel::SRHrecomb3D(int i, int j, int k){

    int pointer = (px)*(py)*(k) + (px)*(j) + (i);
    //http://www.iue.tuwien.ac.at/phd/entner/node11.html
    return (material[pointer].u*material[pointer].v-1)/material[pointer].tau/(material[pointer].u*exp(material[pointer].phi/VT)+material[pointer].v*exp((-1)*material[pointer].phi/VT)+2);
}

double DDmodel::Bern(double dphi, double phi)
{
    if(abs(dphi)<bernXl){
        return exp(phi)*(1.0-dphi/2);
    }
    else{
        return exp(phi)*dphi/(exp(dphi)-1.0);
    }
}

void DDmodel::IdVG2D(){

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
        DDInitialGuess2D();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;

        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=PoissonSolver2D();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=ECSolver2D();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<"\t";

            //hole===============
            errHole=HCSolver2D();
            iter_Hole=DD_loop;
            if(errHole>errMax)
                errMax=errHole;

            output2<<"Hole:" << iter_Hole <<"\t"<<errHole<<"\t";
            output2<<"Max.err:" <<errMax <<endl;

        }while( (iter_Hole!=1 || iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        RhoCalculation2D();
        EfieldCalculation2D();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<".txt";
        name2=name1.str();
        PrintMaterial2D(name2.c_str());

        Jcal2D();

        index++;
        numIter=0;

    }while(volGi+index*volGs<(volGe+0.001));

    output2.close();
    cout << "Simulation Process Finished."<<endl;

}

void DDmodel::IdVD2D(){

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
        DDInitialGuess2D();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;


        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=PoissonSolver2D();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=ECSolver2D();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<"\t";

            //hole===============
            errHole=HCSolver2D();
            iter_Hole=DD_loop;
            if(errHole>errMax)
                errMax=errHole;

            output2<<"Hole:" << iter_Hole <<"\t"<<errHole<<"\t";
            output2<<"Max.err:" <<errMax <<endl;

        }while( (iter_Hole!=1 || iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        RhoCalculation2D();
        EfieldCalculation2D();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<".txt";
        name2=name1.str();
        PrintMaterial2D(name2.c_str());

        Jcal2D();

        index++;
        numIter=0;

    }while(volDi+index*volDs<volDe+0.001);

    output2.close();
    cout << "Simulation Process Finished."<<endl;
}

void DDmodel::IdVG3D(){

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
        DDInitialGuess3D();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;

        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=PoissonSolver3D();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=ECSolver3D();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<"\t";

            //hole===============
            errHole=HCSolver3D();
            iter_Hole=DD_loop;
            if(errHole>errMax)
                errMax=errHole;

            output2<<"Hole:" << iter_Hole <<"\t"<<errHole<<"\t";
            output2<<"Max.err:" <<errMax <<endl;

        }while( (iter_Hole!=1 || iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        RhoCalculation3D();
        EfieldCalculation3D();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<".txt";
        name2=name1.str();
        PrintMaterial3D(name2.c_str());

        Jcal3D();

        index++;
        numIter=0;

    }while(volGi+index*volGs<(volGe+0.001));

    output2.close();
    cout << "Simulation Process Finished."<<endl;

}

void DDmodel::IdVD3D(){

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
        DDInitialGuess3D();

        output2 <<"Vg="<<volG<<"\t"<<"Vd="<<volD<<endl;
        output2 <<"========================================================"<<endl;


        do{
            errMax=0;
            numIter++;
            output2<<numIter<<"\t";

            //poisson============
            errPhi=PoissonSolver3D();
            iter_Phi=DD_loop;
            if(errPhi>errMax)
                errMax=errPhi;

            output2 <<"Poisson:" << iter_Phi <<"\t"<<errPhi<<"\t";

            //electron===========
            errElec=ECSolver3D();
            iter_Elec=DD_loop;
            if(errElec>errMax)
                errMax=errElec;

            output2<<"Electron:" << iter_Elec <<"\t"<<errElec<<"\t";

            //hole===============
            errHole=HCSolver3D();
            iter_Hole=DD_loop;
            if(errHole>errMax)
                errMax=errHole;

            output2<<"Hole:" << iter_Hole <<"\t"<<errHole<<"\t";
            output2<<"Max.err:" <<errMax <<endl;

        }while( (iter_Hole!=1 || iter_Elec!=1 || iter_Phi!=1) && numIter<maxIter);

        output2<<"= = = iteration stop = = ="<<endl<<endl;

        RhoCalculation3D();
        EfieldCalculation3D();

        stringstream name1;
        string name2;

        name1<<"Vg="<<volG<<"_"<<"Vd="<<volD<<".txt";
        name2=name1.str();
        PrintMaterial3D(name2.c_str());

        Jcal3D();

        index++;
        numIter=0;

    }while(volDi+index*volDs<volDe+0.001);

    output2.close();
    cout << "Simulation Process Finished."<<endl;
}


void DDmodel::Jcal2D(){

    switch(StructureFlag){

    case 1:
        Jcal2D_PN();
        break;
    case 2:
        Jcal2D_MOSFET();
        break;
    case 3:
        Jcal2D_ISFET();
        break;
    default:
        cout << "Not appropriate StructureFlag @ Jcal2D."<<endl;
        exit(0);
    }
}

void DDmodel::Jcal2D_PN(void){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0);

    JcalS2D_PN(Current_Sn,Current_Sp);
    JcalD2D_PN(Current_Dn,Current_Dp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific<<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp<<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<endl;
    output1.close();
}

void DDmodel::JcalS2D_PN(double &JSn,double &JSp){

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

        JSn+= (material[pointer].mun+material[pointer].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
             (material[pointer_ip].u-material[pointer].u)*deltay/xstep;

        JSp+= (material[pointer].mup+material[pointer].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
             (material[pointer_ip].v-material[pointer].v)*deltay/xstep;
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::JcalD2D_PN(double &JDn, double &JDp){

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

        JDn+= (material[pointer].mun+material[pointer].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
             (material[pointer_ip].u-material[pointer].u)*deltay/xstep;

        JDp+= (material[pointer].mup+material[pointer].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
             (material[pointer_ip].v-material[pointer].v)*deltay/xstep;
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::Jcal2D_MOSFET(void){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    JcalS2D_MOSFET(Current_Sn,Current_Sp);
    JcalD2D_MOSFET(Current_Dn,Current_Dp);
    JcalB2D_MOSFET(Current_Bn,Current_Bp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific
          <<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
          <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"
          <<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;
    output1.close();
}

void DDmodel::JcalS2D_MOSFET(double &JSn, double &JSp){

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

            JSn+= (material[pointer_jp].mun+material[pointer].mun)/2*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, material[pointer_jp].phi/VT)*
                 (material[pointer_jp].u-material[pointer].u)*deltax/ystep;

            JSp+= (material[pointer_jp].mup+material[pointer].mup)/2*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                 (material[pointer_jp].v-material[pointer].v)*deltax/ystep;
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::JcalD2D_MOSFET(double &JDn, double &JDp){

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

            JDn+= (material[pointer_jp].mun+material[pointer].mun)/2*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, material[pointer_jp].phi/VT)*
                 (material[pointer_jp].u-material[pointer].u)*deltax/ystep;

            JDp+= (material[pointer_jp].mup+material[pointer].mup)/2*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                 (material[pointer_jp].v-material[pointer].v)*deltax/ystep;
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::JcalB2D_MOSFET(double &JBn, double &JBp){

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

        JBn+= (material[pointer_jp].mun+material[pointer].mun)/2*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, material[pointer_jp].phi/VT)*
             (material[pointer_jp].u-material[pointer].u)*deltax/ystep;

        JBp+= (material[pointer].mup+material[pointer].mup)/2*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
             (material[pointer_jp].v-material[pointer].v)*deltax/ystep;
    }

    JBn=JBn*q0*ni_nm*VT*(-1);
    JBp=JBp*q0*ni_nm*VT;
}

void DDmodel::Jcal2D_ISFET(void){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    JcalS2D_ISFET(Current_Sn,Current_Sp);
    JcalD2D_ISFET(Current_Dn,Current_Dp);
    JcalB2D_ISFET(Current_Bn,Current_Bp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific
          <<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
          <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"
          <<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;
    output1.close();
}

void DDmodel::JcalS2D_ISFET(double &JSn, double &JSp){

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

            JSn+= (material[pointer_ip].mun+material[pointer].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                 (material[pointer_ip].u-material[pointer].u)*deltay/xstep;

            JSp+= (material[pointer_ip].mup+material[pointer].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                 (material[pointer_ip].v-material[pointer].v)*deltay/xstep;
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::JcalD2D_ISFET(double &JDn, double &JDp){

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
            JDn+= (material[pointer].mun+material[pointer].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                 (material[pointer_ip].u-material[pointer].u)*deltay/xstep;

            JDp+= (material[pointer].mup+material[pointer].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                 (material[pointer_ip].v-material[pointer].v)*deltay/xstep;
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::JcalB2D_ISFET(double &JBn, double &JBp){

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

        JBn+= (material[pointer_jn].mun+material[pointer].mun)/2*Bern(material[pointer_jn].phi/VT-material[pointer].phi/VT, material[pointer_jn].phi/VT)*
             (material[pointer_jn].u-material[pointer].u)*deltax/ystep;

        JBp+= (material[pointer].mup+material[pointer].mup)/2*Bern(material[pointer_jn].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
             (material[pointer_jn].v-material[pointer].v)*deltax/ystep;
    }

    JBn=JBn*q0*ni_nm*VT*(-1);
    JBp=JBp*q0*ni_nm*VT;

}

void DDmodel::Jcal3D(){

    switch(StructureFlag){

    case 1:
        Jcal3D_PN();
        break;
    case 2:
        Jcal3D_MOSFET();
        break;
    case 3:
        Jcal3D_ISFET();
        break;
    case 4:
        Jcal3D_1NWR();
        break;
    case 5:
        Jcal3D_2NWR();
        break;
    default:
        cout << "Not appropriate StructureFlag @ Jcal3D."<<endl;
        exit(0);
    }
}

double DDmodel::ECSolver3D(){

    DD_loop=0;
    double errEC(0),errEC_max(0);

    do{
        DD_loop++;

        switch(StructureFlag){

        case 1:
            errEC=ECTypeA3D();
            break;
        case 2:
            errEC=ECTypeA3D();
            break;
        case 3:
            errEC=ECTypeB3D();
            break;
        case 4:
            errEC=EC1NWR3D();
            break;
        case 5:
            errEC=EC2NWR3D();
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

double DDmodel::ECTypeA3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<pz-1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = material[pointer].u;

                material[pointer].u=ECInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    ECBC3D();

    return max_val;
}

double DDmodel::HCTypeA3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<pz-1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = material[pointer].v;

                material[pointer].v=HCInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    HCBC3D();

    return max_val;
}

double DDmodel::ECTypeB3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<SubstrateTOP; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = material[pointer].u;

                material[pointer].u=ECInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    ECBC3D();

    return max_val;
}

double DDmodel::HCTypeB3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=1; k<SubstrateTOP; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = material[pointer].v;

                material[pointer].v=HCInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    HCBC3D();

    return max_val;
}

double DDmodel::EC1NWR3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft1+1; i<NWRright1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom1+1; k<NWRtop1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = material[pointer].u;

                material[pointer].u=ECInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    ECBC3D();

    return max_val;

}

double DDmodel::EC2NWR3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft1+1; i<NWRright1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom1+1; k<NWRtop1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double uk = material[pointer].u;

                material[pointer].u=ECInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].u)-log(uk));

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

                double uk = material[pointer].u;

                material[pointer].u=ECInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].u)-log(uk));

                error=error/(abs(VT*log(uk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    ECBC3D_2NWR();

    return max_val;

}

double DDmodel::ECInner3D(int i, int j, int k){

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

    double uf =material[pointer].u;

    double Bip = (material[pointer].mun+material[pointer_ip].mun)/2.0*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT);
    double Bin = (material[pointer].mun+material[pointer_in].mun)/2.0*Bern(material[pointer_in].phi/VT-material[pointer].phi/VT, material[pointer_in].phi/VT);
    double Bjp = (material[pointer].mun+material[pointer_jp].mun)/2.0*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, material[pointer_jp].phi/VT);
    double Bjn = (material[pointer].mun+material[pointer_jn].mun)/2.0*Bern(material[pointer_jn].phi/VT-material[pointer].phi/VT, material[pointer_jn].phi/VT);
    double Bkp = (material[pointer].mun+material[pointer_kp].mun)/2.0*Bern(material[pointer_kp].phi/VT-material[pointer].phi/VT, material[pointer_kp].phi/VT);
    double Bkn = (material[pointer].mun+material[pointer_kn].mun)/2.0*Bern(material[pointer_kn].phi/VT-material[pointer].phi/VT, material[pointer_kn].phi/VT);

    double f=volume/VT*material[pointer].r;

    double df=volume/VT*(pow(material[pointer].v,2)*exp(-material[pointer].phi/VT)+2*material[pointer].v+exp(material[pointer].phi/VT))/material[pointer].tau/
                pow(material[pointer].u*exp(material[pointer].phi/VT)+material[pointer].v*exp(-material[pointer].phi/VT)+2,2);

    uf=((Bip*material[pointer_ip].u/xstep_p+Bin*material[pointer_in].u/xstep_n)*deltay*deltaz
       +(Bjp*material[pointer_jp].u/ystep_p+Bjn*material[pointer_jn].u/ystep_n)*deltax*deltaz
       +(Bkp*material[pointer_kp].u/zstep_p+Bkn*material[pointer_kn].u/zstep_n)*deltax*deltay-f+df*uf)
       /((Bip/xstep_p+Bin/xstep_n)*deltay*deltaz+(Bjp/ystep_p+Bjn/ystep_n)*deltax*deltaz+(Bkp/zstep_p+Bkn/zstep_n)*deltax*deltay+df);

    return uf;
}

double DDmodel::HCSolver3D(){

    DD_loop=0;
    double errHC(0),errHC_max(0);

    do{
        DD_loop++;

        switch(StructureFlag){

        case 1:
            errHC=HCTypeA3D();
            break;
        case 2:
            errHC=HCTypeA3D();
            break;
        case 3:
            errHC=HCTypeB3D();
            break;
        case 4:
            errHC=HC1NWR3D();
            break;
        case 5:
            errHC=HC2NWR3D();
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

double DDmodel::HC1NWR3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft1+1; i<NWRright1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom1+1; k<NWRtop1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = material[pointer].v;

                material[pointer].v=HCInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    HCBC3D();

    return max_val;

}

double DDmodel::HC2NWR3D(){

    double  max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=NWRleft1+1; i<NWRright1; i++) {
        for (int j=1; j<py-1; j++) {
            for (int k=NWRbottom1+1; k<NWRtop1; k++) {

                int pointer = (px)*(py)*(k) + (px)*(j) + (i);

                double vk = material[pointer].v;

                material[pointer].v=HCInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].v)-log(vk));

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

                double vk = material[pointer].v;

                material[pointer].v=HCInner3D(i,j,k);

                material[pointer].r=SRHrecomb3D(i,j,k);

                double error=VT*abs(log(material[pointer].v)-log(vk));

                error=error/(abs(VT*log(vk))+1);

                if(error>max_val)
                    max_val=error;
            }
        }
    }

    HCBC3D_2NWR();

    return max_val;

}

double DDmodel::HCInner3D(int i, int j, int k){

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

    double vf =material[pointer].v;

    double Bip = (material[pointer].mup+material[pointer_ip].mup)/2.0*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);
    double Bin = (material[pointer].mup+material[pointer_in].mup)/2.0*Bern(material[pointer_in].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);
    double Bjp = (material[pointer].mup+material[pointer_jp].mup)/2.0*Bern(material[pointer_jp].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);
    double Bjn = (material[pointer].mup+material[pointer_jn].mup)/2.0*Bern(material[pointer_jn].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);
    double Bkp = (material[pointer].mup+material[pointer_kp].mup)/2.0*Bern(material[pointer_kp].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);
    double Bkn = (material[pointer].mup+material[pointer_kn].mup)/2.0*Bern(material[pointer_kn].phi/VT-material[pointer].phi/VT, (-1)*material[pointer].phi/VT);

    double f=volume/VT*material[pointer].r;

    double df=volume/VT*(pow(material[pointer].u,2)*exp(material[pointer].phi/VT)+2*material[pointer].u+exp(-material[pointer].phi/VT))/material[pointer].tau/
                pow(material[pointer].u*exp(material[pointer].phi/VT)+material[pointer].v*exp(-material[pointer].phi/VT)+2,2);

    vf=((Bip*material[pointer_ip].v/xstep_p+Bin*material[pointer_in].v/xstep_n)*deltay*deltaz
       +(Bjp*material[pointer_jp].v/ystep_p+Bjn*material[pointer_jn].v/ystep_n)*deltax*deltaz
       +(Bkp*material[pointer_kp].v/zstep_p+Bkn*material[pointer_kn].v/zstep_n)*deltax*deltay-f+df*vf)
       /((Bip/xstep_p+Bin/xstep_n)*deltay*deltaz+(Bjp/ystep_p+Bjn/ystep_n)*deltax*deltaz+(Bkp/zstep_p+Bkn/zstep_n)*deltax*deltay+df);

    return vf;
}


void DDmodel::Jcal3D_PN(){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0);

    JcalS3D_PN(Current_Sn,Current_Sp);
    JcalD3D_PN(Current_Dn,Current_Dp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific<<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp<<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<endl;
    output1.close();
}

void DDmodel::JcalS3D_PN(double &JSn,double &JSp){

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

            JSn+= (material[pointer].mun+material[pointer_ip].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                 (material[pointer_ip].u-material[pointer].u)*deltay*deltaz/xstep;

            JSp+= (material[pointer].mup+material[pointer_ip].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer_ip].phi/VT)*
                 (material[pointer_ip].v-material[pointer].v)*deltay*deltaz/xstep;
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::JcalD3D_PN(double &JDn, double &JDp){

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

            JDn+= (material[pointer].mun+material[pointer_ip].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                 (material[pointer_ip].u-material[pointer].u)*deltay*deltaz/xstep;

            JDp+= (material[pointer].mup+material[pointer_ip].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer_ip].phi/VT)*
                 (material[pointer_ip].v-material[pointer].v)*deltay*deltaz/xstep;
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}


void DDmodel::Jcal3D_MOSFET(){

    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    JcalS3D_MOSFET(Current_Sn,Current_Sp);
    JcalD3D_MOSFET(Current_Dn,Current_Dp);
    JcalB3D_MOSFET(Current_Bn,Current_Bp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific
          <<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
          <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"
          <<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;
    output1.close();
}

void DDmodel::JcalS3D_MOSFET(double &JSn, double &JSp){

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

                JSn+= (material[pointer_kp].mun+material[pointer].mun)/2*Bern(material[pointer_kp].phi/VT-material[pointer].phi/VT, material[pointer_kp].phi/VT)*
                     (material[pointer_kp].u-material[pointer].u)*deltax*deltay/zstep;

                JSp+= (material[pointer_kp].mup+material[pointer].mup)/2*Bern(material[pointer_kp].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                     (material[pointer_kp].v-material[pointer].v)*deltax*deltay/zstep;

            }
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::JcalD3D_MOSFET(double &JDn, double &JDp){

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

                JDn+= (material[pointer_kp].mun+material[pointer].mun)/2*Bern(material[pointer_kp].phi/VT-material[pointer].phi/VT, material[pointer_kp].phi/VT)*
                     (material[pointer_kp].u-material[pointer].u)*deltax*deltay/zstep;

                JDp+= (material[pointer_kp].mup+material[pointer].mup)/2*Bern(material[pointer_kp].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                     (material[pointer_kp].v-material[pointer].v)*deltax*deltay/zstep;
            }
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::JcalB3D_MOSFET(double &JBn, double &JBp){

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

            JBn+= (material[pointer_kn].mun+material[pointer].mun)/2*Bern(material[pointer_kn].phi/VT-material[pointer].phi/VT, material[pointer_kn].phi/VT)*
                 (material[pointer_kn].u-material[pointer].u)*deltax*deltay/zstep;

            JBp+= (material[pointer_kn].mup+material[pointer].mup)/2*Bern(material[pointer_kn].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                 (material[pointer_kn].v-material[pointer].v)*deltax*deltay/zstep;
        }
    }

    JBn=JBn*q0*ni_nm*VT*(-1);
    JBp=JBp*q0*ni_nm*VT;
}

void DDmodel::Jcal3D_ISFET(){

    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    JcalS3D_ISFET(Current_Sn,Current_Sp);
    JcalD3D_ISFET(Current_Dn,Current_Dp);
    JcalB3D_ISFET(Current_Bn,Current_Bp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific
          <<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
          <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"
          <<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;
    output1.close();
}

void DDmodel::JcalS3D_ISFET(double &JSn, double &JSp){

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

                JSn+= (material[pointer_ip].mun+material[pointer].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                     (material[pointer_ip].u-material[pointer].u)*deltay*deltaz/xstep;

                JSp+= (material[pointer_ip].mup+material[pointer].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                     (material[pointer_ip].v-material[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::JcalD3D_ISFET(double &JDn, double &JDp){

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

                JDn+= (material[pointer_ip].mun+material[pointer].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                     (material[pointer_ip].u-material[pointer].u)*deltaz*deltay/xstep;

                JDp+= (material[pointer_ip].mup+material[pointer].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                     (material[pointer_ip].v-material[pointer].v)*deltaz*deltay/xstep;
            }
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::JcalB3D_ISFET(double &JBn, double &JBp){

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

            JBn+= (material[pointer_kn].mun+material[pointer].mun)/2*Bern(material[pointer_kn].phi/VT-material[pointer].phi/VT, material[pointer_kn].phi/VT)*
                 (material[pointer_kn].u-material[pointer].u)*deltax*deltay/zstep;

            JBp+= (material[pointer_kn].mup+material[pointer].mup)/2*Bern(material[pointer_kn].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                 (material[pointer_kn].v-material[pointer].v)*deltax*deltay/zstep;
        }
    }

    JBn=JBn*q0*ni_nm*VT*(-1);
    JBp=JBp*q0*ni_nm*VT;
}


void DDmodel::Jcal3D_1NWR(){
    ofstream  output1;
    double Current_Sn(0),Current_Sp(0),Current_Dn(0),Current_Dp(0),Current_Bn(0),Current_Bp(0);

    JcalS3D_NWR1(Current_Sn,Current_Sp);
    JcalD3D_NWR1(Current_Dn,Current_Dp);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<volS<<"\t"<<volG<<"\t"<<volD<<"\t"<<scientific<<Current_Sn<<"\t"<<Current_Sp<<"\t"<<Current_Dn<<"\t"<<Current_Dp
           <<"\t"<<Current_Sn+Current_Sp<<"\t"<<Current_Dn+Current_Dp<<"\t"<<Current_Bn<<"\t"<<Current_Bp<<"\t"<<Current_Bn+Current_Bp<<endl;

    output1.close();
}

void DDmodel::Jcal3D_2NWR(){
    ofstream  output1;
    double Current_Sn1(0),Current_Sp1(0),Current_Dn1(0),Current_Dp1(0);
    double Current_Sn2(0),Current_Sp2(0),Current_Dn2(0),Current_Dp2(0);

    JcalS3D_NWR1(Current_Sn1,Current_Sp1);
    JcalD3D_NWR1(Current_Dn1,Current_Dp1);

    JcalS3D_NWR2(Current_Sn2,Current_Sp2);
    JcalD3D_NWR2(Current_Dn2,Current_Dp2);

    output1.open("current.txt", fstream::out | fstream::app);
    output1.precision(6);
    output1<<(NWRCentery1+NWRCentery2)/2+ShiftDistanceY<<"\t"<<scientific
          <<Current_Sn1<<"\t"<<Current_Sp1<<"\t"<<Current_Dn1<<"\t"<<Current_Dp1<<"\t"<<Current_Sn1+Current_Sp1<<"\t"<<Current_Dn1+Current_Dp1<<"\t"
          <<Current_Sn2<<"\t"<<Current_Sp2<<"\t"<<Current_Dn2<<"\t"<<Current_Dp2<<"\t"<<Current_Sn2+Current_Sp2<<"\t"<<Current_Dn2+Current_Dp2<<"\t"<<endl;

    output1.close();
}

void DDmodel::JcalS3D_NWR1(double &JSn, double &JSp){

    for (int j=NWRleft1; j<NWRright1; j++) {
        for (int k=NWRbottom1; k<NWRtop1; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (0);

            if(material[pointer].Type==1){
                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (1);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (1);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (1);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (1);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;

                JSn+= (material[pointer].mun+material[pointer_ip].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                     (material[pointer_ip].u-material[pointer].u)*deltay*deltaz/xstep;

                JSp+=(material[pointer].mup+material[pointer_ip].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                    (material[pointer_ip].v-material[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::JcalD3D_NWR1(double &JDn, double &JDp){

    for (int j=NWRleft1; j<NWRright1; j++) {
        for (int k=NWRbottom1; k<NWRtop1; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (px-2);
            if(material[pointer].Type==1){
                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (px-1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (px-1);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (px-1);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (px-1);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (px-1);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;

                JDn+= (material[pointer].mun+material[pointer_ip].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                     (material[pointer_ip].u-material[pointer].u)*deltay*deltaz/xstep;

                JDp+=(material[pointer].mup+material[pointer_ip].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                    (material[pointer_ip].v-material[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::JcalS3D_NWR2(double &JSn, double &JSp){

    for (int j=NWRleft1; j<NWRright1; j++) {
        for (int k=NWRbottom1; k<NWRtop1; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (0);

            if(material[pointer].Type==1){
                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (1);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (1);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (1);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (1);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;

                JSn+= (material[pointer].mun+material[pointer_ip].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                     (material[pointer_ip].u-material[pointer].u)*deltay*deltaz/xstep;

                JSp+=(material[pointer].mup+material[pointer_ip].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                    (material[pointer_ip].v-material[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JSn=JSn*q0*ni_nm*VT*(-1);
    JSp=JSp*q0*ni_nm*VT;
}

void DDmodel::JcalD3D_NWR2(double &JDn, double &JDp){

    for (int j=NWRleft1; j<NWRright1; j++) {
        for (int k=NWRbottom1; k<NWRtop1; k++) {

            int pointer = (px)*(py)*(k) + (px)*(j) + (px-2);
            if(material[pointer].Type==1){
                int pointer_ip = (px)*(py)*(k) + (px)*(j) + (px-1);
                int pointer_jp = (px)*(py)*(k) + (px)*(j+1) + (px-1);
                int pointer_jn = (px)*(py)*(k) + (px)*(j-1) + (px-1);
                int pointer_kp = (px)*(py)*(k+1) + (px)*(j) + (px-1);
                int pointer_kn = (px)*(py)*(k-1) + (px)*(j) + (px-1);

                double xstep=abs(mesh[pointer].coordX-mesh[pointer_ip].coordX);
                double deltay=abs(mesh[pointer_jn].coordY-mesh[pointer_jp].coordY)/2;
                double deltaz=abs(mesh[pointer_kn].coordZ-mesh[pointer_kp].coordZ)/2;

                JDn+= (material[pointer].mun+material[pointer_ip].mun)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, material[pointer_ip].phi/VT)*
                     (material[pointer_ip].u-material[pointer].u)*deltay*deltaz/xstep;

                JDp+=(material[pointer].mup+material[pointer_ip].mup)/2*Bern(material[pointer_ip].phi/VT-material[pointer].phi/VT, -material[pointer].phi/VT)*
                    (material[pointer_ip].v-material[pointer].v)*deltay*deltaz/xstep;
            }
        }
    }

    JDn=JDn*q0*ni_nm*VT*(-1);
    JDp=JDp*q0*ni_nm*VT;
}

void DDmodel::FindNWRBC1(){

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

        if(material[pointer].Type==1 && material[pointer_kp].Type==2){
            NWRtop1=k;
            break;
        }
    }

    for(int k=NWRcenterzP1;k>0;k--){

        int pointer = (py)*(px)*(k) + (px)*(NWRcenteryP1) + (px/2);
        int pointer_kn = (py)*(px)*(k-1) + (px)*(NWRcenteryP1) + (px/2);

        if(material[pointer].Type==1 && material[pointer_kn].Type==2){
            NWRbottom1=k;
            break;
        }
    }

    //cerr << "(NWRbottom1,NWRtop1)="<<"("<<NWRbottom1<<","<<NWRtop1<<")"<<endl;

    for(int j=NWRcenteryP1;j<py;j++){

        int pointer = (py)*(px)*(NWRcenterzP1) + (px)*(j) + (px/2);
        int pointer_jp = (py)*(px)*(NWRcenterzP1) + (px)*(j+1) + (px/2);

        if(material[pointer].Type==1 && material[pointer_jp].Type==2){
            NWRright1=j;
            break;
        }
    }

    for(int j=NWRcenteryP1;j<py;j--){

        int pointer = (py)*(px)*(NWRcenterzP1) + (px)*(j) + (px/2);
        int pointer_jn = (py)*(px)*(NWRcenterzP1) + (px)*(j-1) + (px/2);

        if(material[pointer].Type==1 && material[pointer_jn].Type==2){
            NWRleft1=j;
            break;
        }
    }
    //cerr << "(NWRleft1,NWRright1)="<<"("<<NWRleft1<<","<<NWRright1<<")"<<endl;

}

void DDmodel::FindNWRBC2(){

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

        if(material[pointer].Type==1 && material[pointer_kp].Type==2){
            NWRtop2=k;
            break;
        }
    }

    for(int k=NWRcenterzP2;k>0;k--){

        int pointer = (py)*(px)*(k) + (px)*(NWRcenteryP2) + (px/2);
        int pointer_kn = (py)*(px)*(k-1) + (px)*(NWRcenteryP2) + (px/2);

        if(material[pointer].Type==1 && material[pointer_kn].Type==2){
            NWRbottom2=k;
            break;
        }
    }
    //cerr << "(NWRbottom2,NWRtop2)="<<"("<<NWRbottom2<<","<<NWRtop2<<")"<<endl;


    for(int j=NWRcenteryP2;j<py;j++){

        int pointer = (py)*(px)*(NWRcenterzP2) + (px)*(j) + (px/2);
        int pointer_jp = (py)*(px)*(NWRcenterzP2) + (px)*(j+1) + (px/2);

        if(material[pointer].Type==1 && material[pointer_jp].Type==2){
            NWRright2=j;
            break;
        }
    }

    for(int j=NWRcenteryP2;j<py;j--){

        int pointer = (py)*(px)*(NWRcenterzP2) + (px)*(j) + (px/2);
        int pointer_jn = (py)*(px)*(NWRcenterzP2) + (px)*(j-1) + (px/2);

        if(material[pointer].Type==1 && material[pointer_jn].Type==2){
            NWRleft2=j;
            break;
        }
    }
    //cerr << "(NWRleft2,NWRright2)="<<"("<<NWRleft2<<","<<NWRright2<<")"<<endl;
}

void DDmodel::AddDot3D(double DotXCenter, double DotYCenter, int Flag){

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

                    material[pointer].Type=Flag;
                    material[pointer].k=AnalytePermittivity;
                }
            }
        }
    }

    for(int k=0;k<pz;k++){
        for(int j=0;j<py;j++){
            for(int i=0;i<px;i++){

                int pointer= (py)*(px)*(k) + (px)*(j) + (i);

                if(pow(mesh[pointer].coordX-DotXCenter,2)+pow(mesh[pointer].coordY-DotYCenter,2)+pow(mesh[pointer].coordZ-DotZCenter,2)<=AnalyteRadius*AnalyteRadius){
                    material[pointer].rho=AnalyteValence*q0/volume; //q0/volumez
                    //rho = [C/nm3]
                }
            }
        }
    }
}


void DDmodel::AddDotString3D(double Yshift){

    double DotDistance=NWRLength/DotNumber;

    for(int i=0;i<DotNumber;i++){
        AddDot3D(JunctionLength+DotDistance/2+DotDistance*i,(NWRCentery1+NWRCentery2)/2-Yshift,i+1000);
    }
}

