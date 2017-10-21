#include "fvmesh.h"

#include <fstream>

using namespace std;

FVMesh::FVMesh()
{

}

FVMesh::~FVMesh()
{

}

void FVMesh::DDmodelStructurePatameterSet2D(int Struct){

    StructureFlag=Struct;

    fstream output;
    output.open("MeshParameter2D.txt", fstream::out | fstream::trunc);

    cout << "2D Device Simulation."<<endl;

    switch(StructureFlag){

    case 1:
        cout << "Device Structure:PN junction(x direction N:P)." << endl;
        output << "Device Structure:PN junction(x direction N:P)." << endl;

        lx=500;
        ly=100;

        break;

    case 2:
        cout << "Device Structure:MOSFET." << endl;
        output << "Device Structure:MOSFET." << endl;

        lx=500;
        ly=500;

        JunctionLength=lx/5;
        JunctionDepth=ly/10;
        output << "JunctionLength="<<JunctionLength<< endl;
        output << "JunctionDepth="<<JunctionDepth<< endl;
        break;

    case 3:
        cout << "Device Structure:ISFET." << endl;
        output << "Device Structure:ISFET." << endl;

        lx=100;
        SubstrateThickness=100;
        Tox=10;
        ElectrolyteThickness=50;
        JunctionLength=lx/5;
        JunctionDepth=SubstrateThickness/5;
        ly=SubstrateThickness+Tox+ElectrolyteThickness;

        output << "JunctionLength="<<JunctionLength<< endl;
        output << "JunctionDepth="<<JunctionDepth<< endl;
        output << "SubstrateThickness="<<SubstrateThickness<< endl;
        output << "Tox="<<Tox<< endl;

        break;

    default:
        cout << "Undifined Device Structure." << endl;
        output << "Undifined Device Structure." << endl;
        exit(0);
    }
    output.close();
}

void FVMesh::CDStructurePatameterSet2D(){

    cout << "2D CD Simulation."<<endl;

    lx=100*1000;
    ly=100*1000;

}

void FVMesh::DDmodelStructurePatameterSet3D(int Struct){

    StructureFlag=Struct;

    fstream output;
    output.open("MeshParameter3D.txt", fstream::out | fstream::trunc);

    cout << "3D Device Simulation."<<endl;

    switch(StructureFlag){

    case 1:
        cout << "Device Structure:PN junction(x direction N:P)." << endl;
        output << "Device Structure:PN junction(x direction N:P)." << endl;

        lx=500;
        ly=100;
        lz=100;

        break;

    case 2:
        cout << "Device Structure:MOSFET." << endl;
        output << "Device Structure:MOSFET." << endl;

        lx=500;
        ly=100;
        lz=500;
        JunctionLength=lx/5;
        JunctionDepth=lz/10;
        output << "JunctionLength="<<JunctionLength<< endl;
        output << "JunctionDepth="<<JunctionDepth<< endl;
        break;

    case 3:
        cout << "Device Structure:ISFET." << endl;
        output << "Device Structure:ISFET." << endl;

        lx=100;
        SubstrateThickness=100;
        Tox=10;
        ElectrolyteThickness=50;
        JunctionLength=lx/5;
        JunctionDepth=SubstrateThickness/5;
        lz=SubstrateThickness+Tox+ElectrolyteThickness;

        ly=50;

        output << "JunctionLength="<<JunctionLength<< endl;
        output << "SubstrateThickness="<<SubstrateThickness<< endl;
        output << "Tox="<<Tox<< endl;
        break;

    case 4:
        cout << "Device Structure:1NWR." << endl;
        output << "Device Structure:1NWR." << endl;

        NWRradiusy=10/2.0;   //W
        NWRradiusz=10/2.0;   //H
        NWRLength=100;
        JunctionLength=10;
        BOX=100;
        Tox=4;
        SubstrateThickness=100;
        NWRCenterz1=BOX+SubstrateThickness+NWRradiusz;
        ElectrolyteThickness=50;
        JunctionLength=lx/5;
        lz=SubstrateThickness+Tox+ElectrolyteThickness;

        lx=JunctionLength*2+NWRLength;
        ly=60;
        lz=SubstrateThickness+BOX+ElectrolyteThickness;

        NWRCentery1=ly/2;

        output << "JunctionLength="<<JunctionLength<< endl;
        output << "SubstrateThickness="<<SubstrateThickness<< endl;
        output << "BOX="<<BOX<< endl;
        output << "Tox="<<Tox<< endl;
        output << "ElectrolyteThickness="<<ElectrolyteThickness<< endl;
        output << "Tox="<<Tox<< endl;

        output << "NWRLength="<<NWRLength<< " NWRadiusx="<<NWRradiusy<< " NWRadiusz="<<NWRradiusz<< endl;
        output << "NWRCentery1="<<NWRCentery1<< " NWRCenterz1="<<NWRCenterz1<< " JunctionLength="<<JunctionLength<< endl;
        break;

    case 5:
        cout << "Device Structure:2NWR." << endl;
        output << "Device Structure:2NWR." << endl;

        NWRradiusy=10/2.0;   //W
        NWRradiusz=10/2.0;   //H
        NWRLength=100;
        JunctionLength=10;
        BOX=100;
        Tox=10;
        SubstrateThickness=100;
        NWRCenterz1=BOX+SubstrateThickness+NWRradiusz;
        NWRCenterz2=NWRCenterz1;
        ElectrolyteThickness=50;
        lz=SubstrateThickness+Tox+ElectrolyteThickness;

        lx=JunctionLength*2+NWRLength;
        ly=100;
        lz=SubstrateThickness+BOX+ElectrolyteThickness;

        NWRCentery1=ly/4;
        NWRCentery2=ly*3/4;

        output << "JunctionLength="<<JunctionLength<< endl;
        output << "SubstrateThickness="<<SubstrateThickness<< endl;
        output << "BOX="<<BOX<< endl;
        output << "Tox="<<Tox<< endl;
        output << "ElectrolyteThickness="<<ElectrolyteThickness<< endl;

        output << "NWRLength="<<NWRLength<< " NWRadiusx="<<NWRradiusy<< " NWRadiusz="<<NWRradiusz<< endl;
        output << "NWRCentery1="<<NWRCentery1<< " NWRCenterz1="<<NWRCenterz1<< endl;
        output << "NWRCentery2="<<NWRCentery2<< " NWRCenterz2="<<NWRCenterz2<< endl;
        break;

    default:
        cout << "Undifined Device Structure." << endl;
        output << "Undifined Device Structure." << endl;
        exit(0);
    }
    output.close();
}

void FVMesh::DDmodelMeshParameterSet2D(){

    //Dimensions number
    Dimension=2;

    //sensor surface
    //NWRradiusy=lx/20;

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
        meshx[i]=0.5*(i+1);
    }
    /*
    meshx[0]=1e-3;
    meshx[1]=1e-2;
    meshx[2]=1e-3;
    */
    for(int i=0;i<My;i++){
        meshy[i]=0.5*(i+1);
    }

    /*
    meshx[3]=1e-2;
    meshx[4]=1e-3;

    meshy[0]=1;
    meshy[1]=2e-2;
    meshy[2]=1e-3;
    */

    //set initial value(minimum)
    px=py=1;

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
    L=px*py;
}

void FVMesh::CDMeshParameterSet2D(){

    //Dimensions number
    Dimension=2;

    //sensor surface
    //NWRradiusy=lx/20;

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

    //set initial value(minimum)
    px=py=1;

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
    L=px*py;
}

void FVMesh::BlockMeshingMesh2D(){

    mesh = new Mesh [L];

    //assign all coordinate
    MeshInitialize();

    for(int m=0;m<Mx;m++){

        double a= xpin[m];

        for(int i=xb[m];i<xb[m+1]+1;i++){
            for (int j=0;j<py;j++){
                int pointer = (px)*(j) + (i);
                mesh[pointer].coordX=a+(i-xb[m])/meshx[m];
            }
        }
    }

    for(int m=0;m<My;m++){

        double a= ypin[m];

        for (int i=0;i<px;i++){
            for(int j=yb[m];j<yb[m+1]+1;j++){
                int pointer = (px)*(j) + (i);
                mesh[pointer].coordY=a+(j-yb[m])/meshy[m];
            }
        }
    }
}

void FVMesh::DDmodelMeshParameterSet3D(){

    //Dimensions number
    Dimension=3;

    // xyz block numbers.
    Mx=1;
    My=3;
    Mz=3;

    // set xyz pins
    xpin=new double [Mx+1];
    ypin=new double [My+1];
    zpin=new double [Mz+1];

    for(int i=0;i<Mx+1;i++){
        xpin[i]=0+lx/Mx*i;
    }

    /*
    for(int i=0;i<My+1;i++){
        ypin[i]=0+ly/My*i;
    }
    */
    ypin[0]=0;
    ypin[1]=10;
    ypin[2]=90;
    ypin[3]=ly;

    /*
    ypin[0]=0;
    ypin[1]=0;
    ypin[2]=0;
    ypin[3]=0;
    */

    /*
    for(int i=0;i<Mz+1;i++){
        zpin[i]=0+lz/Mz*i;
    }
    */

    zpin[0]=0;
    zpin[1]=190;
    zpin[2]=240;
    zpin[3]=lz;

    // set xyz mesh steps
    meshx=new double [Mx];
    meshy=new double [My];
    meshz=new double [Mz];

    // mesh unit = 1/nm = 1/step
    for(int i=0;i<Mx;i++){
        meshx[i]=1*(i+1);
    }

    /*
    for(int i=0;i<My;i++){
        meshy[i]=1*(i+1);
    }
    */

    meshy[0]=0.4;
    meshy[1]=1;
    meshy[2]=0.4;

    for(int i=0;i<Mz;i++){
        meshz[i]=1*(i+1);
    }

    meshz[0]=0.1;
    meshz[1]=1;
    meshz[2]=0.5;
    /*
    */

    //set initial value(minimum)
    px=py=pz=1;

    // points calculation

    for(int i=0;i<Mx;i++){
        px=px+meshx[i]*(xpin[i+1]-xpin[i]);
    }
    for(int i=0;i<My;i++){
        py=py+meshy[i]*(ypin[i+1]-ypin[i]);
    }
    for(int i=0;i<Mz;i++){
        pz=pz+meshz[i]*(zpin[i+1]-zpin[i]);
    }

    // set xyz  point numbers till each block
    xb=new int [Mx+1];
    yb=new int [My+1];
    zb=new int [Mz+1];

    for(int i=1;i<Mx+1;i++){
        xb[0]=0;
        xb[i]=xb[i-1]+(xpin[i]-xpin[i-1])*meshx[i-1];
    }
    for(int i=1;i<My+1;i++){
        yb[0]=0;
        yb[i]=yb[i-1]+(ypin[i]-ypin[i-1])*meshy[i-1];
    }
    for(int i=1;i<Mz+1;i++){
        zb[0]=0;
        zb[i]=zb[i-1]+(zpin[i]-zpin[i-1])*meshz[i-1];
    }

    L=px*py*pz;
}

void FVMesh::BlockMeshingMesh3D(){

    mesh = new Mesh [L];

    MeshInitialize();

    for(int m=0;m<Mx;m++){

        double a= xpin[m];

        for(int i=xb[m];i<xb[m+1]+1;i++){
            for (int j=0;j<py;j++){
                for (int k=0;k<pz;k++){
                    int pointer =(px)*(py)*(k) + (px)*(j) + (i);
                    mesh[pointer].coordX=a+(i-xb[m])/meshx[m];
                }
            }
        }
    }

    for(int m=0;m<My;m++){

        double a= ypin[m];

        for (int i=0;i<px;i++){
            for(int j=yb[m];j<yb[m+1]+1;j++){
                for (int k=0;k<pz;k++){
                    int pointer =(px)*(py)*(k) + (px)*(j) + (i);
                    mesh[pointer].coordY=a+(j-yb[m])/meshy[m];
                }
            }
        }
    }

    for(int m=0;m<Mz;m++){

        double a= zpin[m];

        for (int i=0;i<px;i++){
            for (int j=0;j<py;j++){
                for(int k=zb[m];k<zb[m+1]+1;k++){
                    int pointer =(px)*(py)*(k) + (px)*(j) + (i);
                    mesh[pointer].coordZ=a+(k-zb[m])/meshz[m];
                }
            }
        }
    }
}

void FVMesh::MeshInitialize(){

    #pragma omp parallel for
    for(int i=0;i<L;i++){
        mesh[i].coordX=0;
        mesh[i].coordY=0;
        mesh[i].coordZ=0;
    }

}

void FVMesh::PrintCoordinate2D(string path){

    //Print Only Coordinate

    path.append(".msh");

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(8);

    output << "X(1)\tY(2)\tZ(3)\t#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            output << mesh[pointer].coordX << '\t' <<  mesh[pointer].coordY  << endl;
        }
    }

    output.close();
}

void FVMesh::PrintCoordinate3D(string path){

    //Print Only Coordinate

    path.append(".msh");

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(8);

    output << "X(1)\tY(2)\tZ(3)\t#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            for (int k=0;k<pz;k++){
                int pointer = (px)*(py)*(k) + (px)*(j) + (i);
                output << mesh[pointer].coordX << '\t' <<  mesh[pointer].coordY  << '\t' << mesh[pointer].coordZ  << endl;
            }
        }
    }
    output.close();
}

void FVMesh::PrintMeshParameter2D(){

    fstream output;
    output.open("MeshParameter2D.txt", fstream::out | fstream::app);
    output << "lx="<<lx<< " ly="<<ly<< endl;
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

    output << "px="<<px<< " py="<<py<< endl<< endl;
    output.close();
}


void FVMesh::PrintMeshParameter3D(){

    fstream output;
    output.open("MeshParameter3D.txt", fstream::out | fstream::app);
    output << "lx="<<lx<< " ly="<<ly<< " lz="<<lz<< endl;
    output << "NWRLength="<<NWRLength<< " NWRadiusx="<<NWRradiusy<< " NWRadiusz="<<NWRradiusz<< endl;
    output << "NWRCentery="<<NWRCentery1<< " NWRCenterz="<<NWRCenterz1<< " JunctionLength="<<JunctionLength<< endl;
    output << "SubstrateThickness="<<SubstrateThickness<< " BOX="<<BOX<< " Tox="<<Tox<< endl;
    output << "Mx="<<Mx<< " My="<<My<< " Mz="<<Mz<< endl;
    for(int i=0;i<Mx+1;i++){
        output <<"xpin["<<i<<"]="<< xpin[i]<<" ";
    }
    output <<endl;
    for(int i=0;i<My+1;i++){
        output <<"ypin["<<i<<"]="<< ypin[i]<<" ";
    }
    output <<endl;
    for(int i=0;i<Mz+1;i++){
        output <<"zpin["<<i<<"]="<< zpin[i]<<" ";
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
    for(int i=0;i<Mz;i++){
        output <<"meshz["<<i<<"]="<< meshz[i]<<" ";
    }
    output <<endl;

    output << "px="<<px<< " py="<<py<< " pz="<<pz<< endl;

    output.close();
}
