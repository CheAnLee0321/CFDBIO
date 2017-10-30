#ifndef FVMESH_H
#define FVMESH_H

#include <iostream>

using namespace std;

struct Mesh {
    double coordX; //coordinate
    double coordY;
    double coordZ;
};

class FVMesh
{
public:
    FVMesh();
    ~FVMesh();


    void FVMesh_DDmodelStructurePatameterSet2D(int Struct);
    void FVMesh_DDmodelStructurePatameterSet3D(int Struct);

    void FVMesh_CDStructurePatameterSet2D();
    void FVMesh_CFDStructurePatameterSet2D();

    //2D
    void FVMesh_DDmodelMeshParameterSet2D();
    void FVMesh_BlockMeshingMesh2D();
    void FVMesh_CDMeshParameterSet2D();
    void FVMesh_CFDMeshParameterSet2D();

    //3D
    void FVMesh_DDmodelMeshParameterSet3D();
    void FVMesh_BlockMeshingMesh3D();

    //Tool function
    void FVMesh_PrintCoordinate2D(string path);
    void FVMesh_PrintCoordinate3D(string path);
    void FVMesh_PrintMeshParameter2D();
    void FVMesh_PrintMeshParameter3D();

private:
    //Tool function
    void FVMesh_MeshInitialize();

protected:

    /*
    StructureFlag=1 PN
    StructureFlag=2 MOSFET
    StructureFlag=3 ISFET
    StructureFlag=4 1NWR in 3D
    StructureFlag=5 2NWR in 3D
    StructureFlag=0 undefined structure
    */

    int StructureFlag=0;
    int Dimension=0;

    //*********************************
    //Structure Build
    //*********************************
    Mesh *mesh;

    double lx=0,ly=0,lz=0,*xpin,*ypin,*zpin,*meshx,*meshy,*meshz;
    int Mx=0,My=0,Mz=0,px=0,py=0,pz=0,*xb,*yb,*zb,L;
    double NWRradiusy=0,NWRradiusz=0,NWRLength=0,NWRCentery1=0,NWRCenterz1=0,NWRCentery2=0,NWRCenterz2=0;

    double JunctionDepth=0, JunctionLength=0, Tox=0;
    double SubstrateThickness=0, BOX=0, ChannelLength=0, ElectrolyteThickness=0;

};

#endif // FVMESH_H
