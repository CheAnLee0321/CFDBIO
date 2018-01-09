//#include "fvmesh.h"
#include "ddmodel.h"
#include "montecarlo.h"
//#include "cdmodel.h"

//#include "cdmodel.h"
#include <iostream>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int main()
{


    //Apply 2D CFD to 3D particle tracing
    MonteCarlo *test=new MonteCarlo();
    test->FVMesh_CFDStructurePatameterSet2D();
    test->FVMesh_CFDMeshParameterSet2D();
    test->FVMesh_BlockMeshingMesh2D();
    test->FVMesh_PrintMeshParameter2D();
    test->CFD_NSParameter();
    test->CFD_NewAndInitialize();
    cout << "Reading Field Data." << endl;
    //test->CFD_ReadVelocity2D("SIMPLER_2D.txt");
    //test->CFD_PrintVelocity2D("Read.txt");

    cout << "Start Monte-Carlo simulation." <<endl;
    test->MC_ParticleTracingParameter3D();
    test->MC_ParticleTracingNew();
    test->MC_ParticleTracingInitialize3D();
    test->MC_ParticleTracingSimulation3D();
    /*
    */

    /*
    //2D CFD template
    CFD *test = new CFD();
    test->FVMesh_CFDStructurePatameterSet2D();
    test->FVMesh_CFDMeshParameterSet2D();
    test->FVMesh_BlockMeshingMesh2D();
    test->FVMesh_PrintMeshParameter2D();
    test->CFD_NSParameter();
    //test->CFD_PrintNSParameter();
    test->CFD_NewAndInitialize();
    //test->CFD_ACPoissonInitialGuess2D();
    //test->CFD_PrintMaterialComplex2D("Initial.txt");
    //test->CFD_PoissonSolverComplex2D();
    //test->CFD_PrintPotential2D("Potential2D.txt");
    test->CFD_ReadPotential2D("Potential2D.txt");
    test->CFD_NSInitialGuessComplex2D();
    test->CFD_VelocityCalculation2D();
    test->CFD_MaxSpeed2D();
    test->CFD_PrintVelocity2D("InitialVelocity.txt");
    test->CFD_SIMPLER2D();
    */

    /*
    //3D CFD template
    CFD *test = new CFD();
    test->FVMesh_CFDStructurePatameterSet3D();
    test->FVMesh_CFDMeshParameterSet3D();
    test->FVMesh_BlockMeshingMesh3D();
    test->FVMesh_PrintMeshParameter3D();
    test->CFD_NSParameter();
    test->CFD_NewAndInitialize();
    test->CFD_ACPoissonInitialGuess3D();
    test->CFD_PrintMaterialComplex3D("Initial.txt");
    test->CFD_PoissonSolverComplex3D();
    test->CFD_PrintPotential3D("Potential.txt");
    test->CFD_NSInitialGuessComplex3D();
    test->CFD_VelocityCalculation3D();
    test->CFD_MaxSpeed();
    test->CFD_PrintVelocity3D("InitialVelocity.txt");
    test->CFD_SIMPLER3D();
    */


    //2D CFD template
    /*
    CFD *test = new CFD();
    test->FVMesh_CFDStructurePatameterSet2D();
    test->FVMesh_CFDMeshParameterSet2D();
    test->FVMesh_BlockMeshingMesh2D();
    test->FVMesh_PrintMeshParameter2D();
    test->CFD_NSParameter();
    test->CFD_NewInitialize();
    test->CFD_ACPoissonInitialGuess2D();
    test->CFD_PrintMaterialComplex2D("I.txt");
    test->CFD_PoissonSolverComplex2D();
    test->CFD_ReadMaterialComplex2D("P.txt");
    test->CFD_PrintMaterialComplex2D("PP.txt");
    test->CFD_NSInitialGuessComplex2D();
    test->CFD_VelocityCalculation2D();
    test->CFD_PrintMaterialComplex2D("CFDI.txt");
    test->CFD_SIMPLER2D();
    */

    /*
    //Monte-Carlo LOD template
    MonteCarlo *test=new MonteCarlo();
    test->MC_ParticleTracingParameter3D();
    test->MC_ParticleTracingNew();
    test->MC_ParticleTracingInitialize3D();
    test->MC_DirectGenerateOnSurface();
    test->MC_Distribution();
    */

    /*
    //2D Convection Diffusion model template
    CDmodel *test=new CDmodel();
    test->FVMesh_CDStructurePatameterSet2D();
    test->FVMesh_CDMeshParameterSet2D();
    test->FVMesh_BlockMeshingMesh2D();
    test->FVMesh_PrintCoordinate2D("test");
    test->CD_Parameter();
    test->CD_AddReceptors2D();
    test->CD_InitialGuess2D();
    test->CD_Solver2D();
    test->CD_PrintMaterial2D("E.txt");
    */

    //MC template
    /*
    MonteCarlo *test=new MonteCarlo();
    test->MC_ParticleTracingParameter3D();
    test->MC_ParticleTracingNew();
    test->MC_ParticleTracingInitialize3D();
    //test->MC_DirectGenerateOnSurface();
    test->MC_ParticleTracingSimulation3D();
    */

    /*
    //2D DD model template
    DDmodel *test2DPN=new DDmodel;
    test2DPN->FVMesh_DDmodelStructurePatameterSet2D(3);
    test2DPN->FVMesh_DDmodelMeshParameterSet2D();
    test2DPN->FVMesh_BlockMeshingMesh2D();
    test2DPN->FVMesh_PrintCoordinate2D("2D");
    test2DPN->FVMesh_PrintMeshParameter2D();
    test2DPN->DD_ParameterSet();
    test2DPN->DD_NewAndInitialize();
    test2DPN->DD_InitialGuess2D();
    //test2DPN->PrintMaterial2D("ISFET.txt");
    //test2DPN->IdVD2D();
    test2DPN->DD_IdVG2D();
    */


    cout << "Simulation Success." << endl;
    return 0;
}
