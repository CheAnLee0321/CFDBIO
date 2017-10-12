#include "ddmodel.h"
#include "fvmesh.h"
#include "montecarlo.h"

#include <iostream>

using namespace std;

int main()
{

    MonteCarlo *test=new MonteCarlo();
    test->ParticleTracingParameter3D();
    test->ParticleTracingNew();
    test->ParticleTracingInitialize3D();
    //test->DirectGenerateOnSurface();
    test->ParticleTracingSimulation3D();
    /*
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
