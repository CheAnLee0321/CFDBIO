//#include "fvmesh.h"
//#include "ddmodel.h"
#include "montecarlo.h"

//#include "cdmodel.h"
#include <iostream>

using namespace std;

int main()
{

    CFD *test = new CFD();
    test->FVMesh_CFDStructurePatameterSet2D();
    test->FVMesh_CFDMeshParameterSet2D();
    test->FVMesh_BlockMeshingMesh2D();
    test->CFD_NSParameter();
    test->CFD_ACPoissonInitialGuess2D();
    test->CFD_PrintMaterialComplex2D("I.dat");
    test->CFD_PoissonSolverComplex2D();
    //test->CFD_NSInitialGuessComplex2D();
    //test->CFD_VelocityCalculation2D();
    //test->PrintCoordinate2D("C.msh");
    test->CFD_PrintMaterialComplex2D("M.dat");


    //MonteCarlo *test=new MonteCarlo();
    //test->ParticleTracingParameter3D();
    //test->ParticleTracingNew();
    //test->ParticleTracingInitialize3D();
    //test->DirectGenerateOnSurface();
    //test->Distribution();


    /*
    CDmodel *test=new CDmodel();
    test->CDStructurePatameterSet2D();
    test->CDMeshParameterSet2D();
    test->BlockMeshingMesh2D();
    test->PrintCoordinate2D("test");
    test->CDParameter();
    test->CDAddReceptors2D();
    test->CDInitialGuess2D();
    test->CDSolver2D();
    test->PrintMaterial2D("E.txt");
    */

    /*
    MonteCarlo *test=new MonteCarlo();
    test->ParticleTracingParameter3D();
    test->ParticleTracingNew();
    test->ParticleTracingInitialize3D();
    //test->DirectGenerateOnSurface();
    test->ParticleTracingSimulation3D();
    */

    /*
    DDmodel *test2DPN=new DDmodel;
    test2DPN->DDmodelStructurePatameterSet2D(3);
    test2DPN->DDmodelMeshParameterSet2D();
    test2DPN->BlockMeshingMesh2D();
    test2DPN->PrintCoordinate2D("2D");
    test2DPN->PrintMeshParameter2D();

    test2DPN->DDmodelParameterSet();
    test2DPN->DDInitialGuess2D();
    //test2DPN->PrintMaterial2D("ISFET.txt");
    //test2DPN->IdVD2D();
    test2DPN->IdVG2D();
    */


    cout << "Simulation Success." << endl;
    return 0;
}
