#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSimpleTimer.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzlog.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <DarcyFlow/TPZDarcyFlow.h>
#include <TPZGenGrid2D.h>
#include <TPZGenGrid3D.h>
#include <TPZLinearAnalysis.h>
#include <iostream>

using namespace std;

enum EnumMatids {
  EMatId = 1, //-------------------------NAO ENTENDI ESSES VALORES-----------------------------
  EBottom = 2,
  ERight = 3,
  ETop = 4,
  ELeft = 5
};

TPZGeoMesh *createMeshWithGenGrid(const TPZVec<int> &nelDiv, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX) {
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  TPZGenGrid2D generator(nelDiv, minX, maxX);
  generator.SetElementType(MMeshType::EQuadrilateral);
  generator.Read(gmesh, EMatId);
  generator.SetBC(gmesh, 4, EBottom);
  generator.SetBC(gmesh, 5, ERight);
  generator.SetBC(gmesh, 6, ETop);
  generator.SetBC(gmesh, 7, ELeft);

  return gmesh;
}

TPZGeoMesh *createRectangularGmesh() {
  TPZGeoMesh *gmesh = new TPZGeoMesh;

  TPZManVector<REAL, 3> coord(3, 0.); // cria um vetor de 3 posicoes inicializadas com 0.

  TPZGeoNode nod0(0, coord, *gmesh); //(0,0)
  coord[0] = 1.;
  TPZGeoNode nod1(1, coord, *gmesh); //(1,0)
  coord[1] = 1.;
  TPZGeoNode nod2(2, coord, *gmesh); //(1,1)
  coord[0] = 0.;
  TPZGeoNode nod3(3, coord, *gmesh); //(0,1)
  coord[0] = 2.;
  coord[1] = 0.;
  TPZGeoNode nod4(4, coord, *gmesh); //(2,0)
  coord[1] = 1.;
  TPZGeoNode nod5(5, coord, *gmesh); //(2,1)

  gmesh->NodeVec().Resize(6);
  gmesh->NodeVec()[0] = nod0;
  gmesh->NodeVec()[1] = nod1;
  gmesh->NodeVec()[2] = nod2;
  gmesh->NodeVec()[3] = nod3;
  gmesh->NodeVec()[4] = nod4;
  gmesh->NodeVec()[5] = nod5;

  // ------ First quadrilateral element ------
  TPZManVector<int64_t, 4> nodeIndexes(4);
  nodeIndexes[0] = 0;
  nodeIndexes[1] = 1;
  nodeIndexes[2] = 2;
  nodeIndexes[3] = 3;

  int64_t index = 0;
  TPZGeoEl *gel = gmesh->CreateGeoElement(EQuadrilateral, nodeIndexes, EMatId, index);

  // ------ Second quadrilateral element ------
  nodeIndexes[0] = 1;
  nodeIndexes[1] = 4;
  nodeIndexes[2] = 5;
  nodeIndexes[3] = 2;
  index = 1;

  TPZGeoEl *gel2 = gmesh->CreateGeoElement(EQuadrilateral, nodeIndexes, EMatId, index);

  // ------ Creating BCs ------
  gel->CreateBCGeoEl(4, EBottom); // side 4 is the bottom side of the element
  gel2->CreateBCGeoEl(4, EBottom);
  gel2->CreateBCGeoEl(5, ERight);
  gel2->CreateBCGeoEl(6, ETop);
  gel->CreateBCGeoEl(6, ETop);
  gel->CreateBCGeoEl(7, ELeft);

  gmesh->SetDimension(2);
  gmesh->BuildConnectivity();

  return gmesh;
}
// Neighbours for side   0 : 7/1 2/0: o lado 0 é vizinho do lado 1 do elemento 7 e do lado 0 do elemento 2 ->Ver desenho do giovani
// Neighbours for side   1 : 3/0 1/0 2/1  o lado 1 é vizinho do elemento 3 pelo lado 0, do elemento 1 pelo lado 0 e do elemento 2 pelo lado 1

TPZCompMesh *createCompMesh(TPZGeoMesh *gmesh) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(1);
  cmesh->SetAllCreateFunctionsContinuous();

  // Add materials (weak formulation)

  TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
  mat->SetConstantPermeability(1.0);
  cmesh->InsertMaterialObject(mat);

  // Add boundary conditions

  int diritype = 0, neumanntype = 1, mixedtype = 2;
  TPZManVector<REAL, 1> val2(1, 3.); // Tudo o que vai para o vetor de carga
  TPZFMatrix<REAL> val1(1, 1, 0.);   // Tudo o que vai para matriz de rigidez

  // Creating boundary conditions
  val2[0] = 0.;
  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.;
  bcond = mat->CreateBC(mat, ERight, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 0.;
  bcond = mat->CreateBC(mat, ETop, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 3.;
  bcond = mat->CreateBC(mat, ELeft, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  // Set up the computational mesh
  cmesh->AutoBuild();

  return cmesh;
}

int main(int argc, char *const argv[]) {
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  //---------------------------Create geometric mesh---------------------------
  const int nthreads = 24;
  const bool isUseGenGrid = true;
  TPZGeoMesh *gmesh = nullptr;
  if (isUseGenGrid) {
    const int neldiv = 20;
    gmesh = createMeshWithGenGrid({neldiv, neldiv}, {0., 0.}, {2., 1.});
  } else {
    gmesh = createRectangularGmesh();
  }

  gmesh->Print(std::cout);
  std::ofstream out("geomesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  //------------------------Create computational mesh-------------------------

  TPZCompMesh *cmesh = createCompMesh(gmesh);
  cmesh->Print(std::cout);

  //-------------------------Create analysis object--------------------------

  TPZLinearAnalysis an(cmesh);
  TPZSSpStructMatrix<STATE> matsp(cmesh);
  // TPZSkylineStructMatrix<STATE> matsp(cmesh);
  matsp.SetNumThreads(nthreads);
  an.SetStructuralMatrix(matsp);
  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky);
  an.SetSolver(step);
  an.Run();
  an.Solution().Print("Solution");

  // Post-processing
  const std::string plotfile("PostProcess");
  constexpr int vtkRes(0);

  TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.Do();

  return 0;
}
