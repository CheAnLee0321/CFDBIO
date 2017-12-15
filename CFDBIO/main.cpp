//#include "fvmesh.h"
//#include "ddmodel.h"
#include "montecarlo.h"

//#include "cdmodel.h"
#include <iostream>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int main()
{

    //MonteCarlo *test = new MonteCarlo;
    //test->MC_ParticleTracingParameter3D();

    CFD *test = new CFD();
    test->FVMesh_CFDStructurePatameterSet3D();
    test->FVMesh_CFDMeshParameterSet3D();
    test->FVMesh_BlockMeshingMesh3D();
    test->FVMesh_PrintMeshParameter3D();
    test->CFD_NSParameter();
    test->CFD_ACPoissonInitialGuess3D();
    test->CFD_PrintMaterialComplex3D("Initial.dat");
    test->CFD_PoissonSolverComplex3D();
    //test->CFD_PrintMaterialComplex3D("AfterPoisson.dat");
    //test->CFD_ReadMaterialComplex3D("AfterPoisson.dat");
    test->CFD_PrintPotential3D("Potential.dat");
    //test->CFD_ReadPotential3D("Potential.dat");
    test->CFD_NSInitialGuessComplex3D();
    test->CFD_VelocityCalculation3D();
    test->CFD_MaxSpeed();
    test->CFD_PrintVelocity3D("InitialVelocity.dat");
    //test->CFD_ReadVelocity3D("InitialVelocity.dat");
    test->CFD_SIMPLER3D();


    /*
    CFD *test = new CFD();
    test->FVMesh_CFDStructurePatameterSet2D();
    test->FVMesh_CFDMeshParameterSet2D();
    test->FVMesh_BlockMeshingMesh2D();
    test->FVMesh_PrintMeshParameter2D();
    test->CFD_NSParameter();
    test->CFD_ACPoissonInitialGuess2D();
    test->CFD_PrintMaterialComplex2D("I.dat");
    test->CFD_PoissonSolverComplex2D();
    test->CFD_ReadMaterialComplex2D("P.dat");
    test->CFD_PrintMaterialComplex2D("PP.dat");
    test->CFD_NSInitialGuessComplex2D();
    test->CFD_VelocityCalculation2D();
    test->CFD_PrintMaterialComplex2D("CFDI.dat");
    test->CFD_SIMPLER2D();
    */

    /*
    //MonteCarlo *test=new MonteCarlo();
    //test->ParticleTracingParameter3D();
    //test->ParticleTracingNew();
    //test->ParticleTracingInitialize3D();
    //test->DirectGenerateOnSurface();
    //test->Distribution();
    */

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
    system("pause");
    return 0;
}
