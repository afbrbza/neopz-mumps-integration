#include "OutputCapture.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSSpStructMatrixMumps.h"
#include "TPZSYSMPPardiso.h"
#include "TPZSimpleTimer.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzfmatrix.h"
#include "pzlog.h"
#include "pzskylmat.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <DarcyFlow/TPZDarcyFlow.h>
#include <TPZGenGrid2D.h>
#include <TPZGenGrid3D.h>
#include <TPZLinearAnalysis.h>
#include <TPZSYSMPMumps.h>
#include <iostream>

// Simple constant forcing function for the domain
void ForcingFunction(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {
  f.resize(1);
  f[0] = 2.0; // Constant source term
}

enum EnumMatids {
  EMatId = 1,
  EBottom = 2,
  ERight = 3,
  ETop = 4,
  ELeft = 5
};

TPZGeoMesh *createMeshWithGenGrid(const TPZVec<int> &nelDiv, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX) { // ----------------O QUE É O GENGRID?----------------
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

TPZGeoMesh *createMesh3D(const TPZVec<int> &nelDiv) {
  TPZVec<REAL> minX(3, 0.);
  TPZVec<REAL> maxX(3, 1.);
  MMeshType elType = MMeshType::ETetrahedral;
  TPZGenGrid3D gen3d(minX, maxX, nelDiv, elType);
  gen3d.BuildVolumetricElements(1);
  TPZGeoMesh *gmesh = gen3d.BuildBoundaryElements(-1, -3, -1, -2, -1, -1);
  return gmesh;
}

TPZCompMesh *createCompMesh3D(TPZGeoMesh *gmesh, const int pord) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(pord);
  cmesh->SetAllCreateFunctionsContinuous();

  // Add materials (weak formulation)

  TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
  mat->SetConstantPermeability(1.0);
  // mat->SetForcingFunction(ForcingFunction, 2); // Add volume source term
  cmesh->InsertMaterialObject(mat);

  // Add boundary conditions

  int diritype = 0, neumanntype = 1, mixedtype = 2;
  TPZManVector<REAL, 1> val2(1, 3.); // Tudo o que vai para o vetor de carga
  TPZFMatrix<REAL> val1(1, 1, 0.);   // Tudo o que vai para matriz de rigidez

  // Creating boundary conditions
  val2[0] = 0.; // Non-zero Neumann on bottom
  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, -1, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.; // Non-zero Neumann on right
  bcond = mat->CreateBC(mat, -2, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 3; // Non-zero Neumann on left
  bcond = mat->CreateBC(mat, -3, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  // Set up the computational mesh
  cmesh->AutoBuild();

  return cmesh;
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
  int iddiri = 2;
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
// Neighbours for side   0 : 7/1 2/0: o lado 0 é vizinho do lado 1 do elemento 7 e do lado 0 do elemento 2 ->Ver desenho do giovani
// Neighbours for side   1 : 3/0 1/0 2/1  o lado 1 é vizinho do elemento 3 pelo lado 0, do elemento 1 pelo lado 0 e do elemento 2 pelo lado 1

TPZCompMesh *createCompMesh(TPZGeoMesh *gmesh, const int pord) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(pord);
  cmesh->SetAllCreateFunctionsContinuous();

  // Add materials (weak formulation)

  TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
  mat->SetConstantPermeability(1.0);
  // mat->SetForcingFunction(ForcingFunction, 2); // Add volume source term
  cmesh->InsertMaterialObject(mat);

  // Add boundary conditions

  int diritype = 0, neumanntype = 1, mixedtype = 2;
  TPZManVector<REAL, 1> val2(1, 3.); // Tudo o que vai para o vetor de carga
  TPZFMatrix<REAL> val1(1, 1, 0.);   // Tudo o que vai para matriz de rigidez

  // Creating boundary conditions
  val2[0] = 0.; // Non-zero Neumann on bottom
  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.; // Non-zero Neumann on right
  bcond = mat->CreateBC(mat, ERight, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 0; // Non-zero Neumann on top
  bcond = mat->CreateBC(mat, ETop, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 3; // Non-zero Neumann on left
  bcond = mat->CreateBC(mat, ELeft, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);

  // Set up the computational mesh
  cmesh->AutoBuild();

  return cmesh;
}

void runMumps(TPZAutoPointer<TPZCompMesh> cmesh, int neldiv, int pOrder, int nthreads) {
  //-------------------------Create analysis object--------------------------
  TPZLinearAnalysis an(cmesh, RenumType::ENone);
  TPZSSpStructMatrixMumps<STATE> matsp(cmesh);
  // TPZSkylineStructMatrix<STATE> matsp(cmesh);
  matsp.SetNumThreads(nthreads);
  an.SetStructuralMatrix(matsp);

  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky);

  an.SetSolver(step);

  an.Assemble();

  auto *smp = dynamic_cast<TPZSYsmpMatrix<STATE> *>(an.MatrixSolver<STATE>().Matrix().operator->());
  if (!smp) {
    DebugStop();
    return;
  }
  const auto &IA = smp->IA();
  const auto &JA = smp->JA();
  const auto &A = smp->A();
  const int64_t nrows = smp->Rows();
  const int64_t nnz = IA[nrows];

  // auto sparseMatrix = dynamic_cast<TPZSYsmpMatrix<STATE> *>(stiff->Create());
  // auto mCast = mSolverCast.Matrix();
  // auto &mSolverCast = an.MatrixSolver<STATE>();
  // auto mCast = mSolverCast.Matrix();
  // auto &pardisoControl = dynamic_cast<TPZSYsmpMatrixPardiso<STATE> *>(mCast.operator->())->GetPardisoControl();
  // pardisoControl.SetMessageLevel(1);

  auto &mSolverCast = an.MatrixSolver<STATE>();
  auto mCast = mSolverCast.Matrix();
  mCast->SetDefPositive(true); // Say to MUMPS that the matrix is positive definite (if it is) to enable optimizations. If not sure, leave as false.
  auto &solverControl = dynamic_cast<TPZSYsmpMatrixMumps<STATE> *>(mCast.operator->())->GetMumpsControl();
  solverControl.SetMessageLevel(2);

  std::cout << "Number of non-zeros: " << nnz << "\n";

  SolverMetrics metrics;
  std::string capturedOutput;
  double neopzSolveTime = 0.0;
  {
    OutputCapture capture;

    TPZSimpleTimer tPZSolve;
    an.Solve();
    neopzSolveTime = tPZSolve.ReturnTimeDouble() / 1000.0; // convert to seconds

    capturedOutput = capture.getOutput();
    metrics = capture.parseMumpsOutput(capturedOutput);
  }

  // std::cout << capturedOutput << std::endl;

  std::cout << "MUMPS Solver Metrics:" << std::endl
            << "nelDiv: " << neldiv << std::endl
            << "pOrder: " << pOrder << std::endl
            << "Number of equations: " << cmesh->NEquations() << std::endl
            << "Estimated number of non-zeros: " << nnz << std::endl;
  if (metrics.numThreads)
    std::cout << "nThreads: " << metrics.numThreads << std::endl;
  if (metrics.orderingBasedOn.size() > 0)
    std::cout << "Ordering based on: " << metrics.orderingBasedOn << std::endl;
  if (metrics.factorizationTime)
    std::cout << "factorizationTime: " << metrics.factorizationTime << " s" << std::endl;
  if (metrics.solveTime)
    std::cout << "solveTime: " << std::fixed << std::setprecision(7) << metrics.solveTime << " s" << std::endl;
  std::cout << "neopzSolveTime: " << std::fixed << std::setprecision(7) << neopzSolveTime << " s" << std::endl;

  an.Solution().Print("Solution from MUMPS");

  // ---------------------------------------

  // Post-processing
  const std::string plotfile("PostProcessMumps");
  constexpr int vtkRes(0);

  TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.Do();
}

void runPardiso(TPZAutoPointer<TPZCompMesh> cmesh, int neldiv, int pOrder, int nthreads) {
  //-------------------------Create analysis object--------------------------
  TPZLinearAnalysis an(cmesh, RenumType::ENone);
  TPZSSpStructMatrix<STATE> matsp(cmesh);
  // TPZSkylineStructMatrix<STATE> matsp(cmesh);
  matsp.SetNumThreads(nthreads);
  an.SetStructuralMatrix(matsp);

  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky);

  an.SetSolver(step);

  an.Assemble();

  auto *smp = dynamic_cast<TPZSYsmpMatrix<STATE> *>(an.MatrixSolver<STATE>().Matrix().operator->());
  if (!smp) {
    DebugStop();
    return;
  }
  const auto &IA = smp->IA();
  const auto &JA = smp->JA();
  const auto &A = smp->A();
  const int64_t nrows = smp->Rows();
  const int64_t nnz = IA[nrows];

  auto &mSolverCast = an.MatrixSolver<STATE>();
  auto mCast = mSolverCast.Matrix();
  mCast->SetDefPositive(true); // Say to Pardiso that the matrix is positive definite (if it is) to enable optimizations. If not sure, leave as false.
  auto &solverControl = dynamic_cast<TPZSYsmpMatrixPardiso<STATE> *>(mCast.operator->())->GetPardisoControl();
  solverControl.SetMessageLevel(1);

  std::cout << "Number of non-zeros: " << nnz << "\n";

  SolverMetrics metrics;
  std::string capturedOutput;
  double neopzSolveTime = 0.0;
  {
    OutputCapture capture;

    TPZSimpleTimer tPZSolve;
    an.Solve();
    neopzSolveTime = tPZSolve.ReturnTimeDouble() / 1000.0; // convert to seconds

    capturedOutput = capture.getOutput();
    metrics = capture.parsePardisoOutput(capturedOutput);
  }

  // std::cout << capturedOutput << std::endl;

  std::cout << "PARDISO Solver Metrics:" << std::endl
            << "nelDiv: " << neldiv << std::endl
            << "pOrder: " << pOrder << std::endl
            << "Number of equations: " << cmesh->NEquations() << std::endl
            << "Estimated number of non-zeros: " << nnz << std::endl;
  if (metrics.numThreads)
    std::cout << "nThreads: " << metrics.numThreads << std::endl;
  if (metrics.orderingBasedOn.size() > 0)
    std::cout << "Ordering based on: " << metrics.orderingBasedOn << std::endl;
  if (metrics.factorizationTime)
    std::cout << "factorizationTime: " << metrics.factorizationTime << " s" << std::endl;
  if (metrics.solveTime)
    std::cout << "solveTime: " << std::fixed << std::setprecision(7) << metrics.solveTime << " s" << std::endl;
  std::cout << "neopzSolveTime: " << std::fixed << std::setprecision(7) << neopzSolveTime << " s" << std::endl;

  an.Solution().Print("Solution from Pardiso");

  // ---------------------------------------

  // Post-processing
  const std::string plotfile("PostProcessPardiso");
  constexpr int vtkRes(0);

  TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.Do();
}

int main(int argc, char *const argv[]) {
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  //---------------------------Create geometric mesh---------------------------
  const int nthreads = 16;
  const bool isUseGenGrid = true;
  TPZGeoMesh *gmesh = nullptr;
  int neldiv = 0, n;
  const bool is3d = true;
  if (isUseGenGrid) {
    if (argc > 1) {
      neldiv = atoi(argv[1]);
    } else {
      std::cout << "Enter number of divisions: ";
      std::cin >> n;
      // neldiv = (std::sqrt(n) - 1);
      neldiv = n;
    }
    std::cout << "Number of divisions in each direction: " << neldiv << std::endl;

    // std::cout<< neldiv;
    std::cout << "Total number of nodes: " << ((neldiv + 1) * (neldiv + 1)) << std::endl;
    if (neldiv <= 0) {
      std::cerr << "Invalid input, using 1 division." << std::endl;
      neldiv = 1;
    }
    if (is3d) {
      // neldiv = 5;
      gmesh = createMesh3D({neldiv, neldiv, neldiv});
    } else
      gmesh = createMeshWithGenGrid({neldiv, neldiv}, {0., 0.}, {1., 1.});
  } else {
    gmesh = createRectangularGmesh();
  }

  // gmesh->Print(std::cout);
  // std::ofstream out("geomesh.vtk");
  // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  //------------------------Create computational mesh-------------------------
  const int pOrder = argc > 2 ? atoi(argv[2]) : 1;
  const int pord = pOrder; // Polynomial order
  TPZAutoPointer<TPZCompMesh> cmesh = nullptr;
  if (is3d) {
    cmesh = createCompMesh3D(gmesh, pord);
  } else {
    cmesh = createCompMesh(gmesh, pord);
  }

  std::cout << "Number of equations: " << cmesh->NEquations() << "\n";

  // TPZAutoPointer<TPZCompMesh> cmesh = createCompMesh(gmesh, pord);
  // cmesh->Print(std::cout);

  #ifdef PZ_USING_MKL
  {
    TPZAutoPointer<TPZCompMesh> cmeshPardiso = new TPZCompMesh(*cmesh);
    runPardiso(cmeshPardiso, neldiv, pord, nthreads);
  }
  #endif
  
  #ifdef PZ_USING_MUMPS
  {
    TPZAutoPointer<TPZCompMesh> cmeshMumps = new TPZCompMesh(*cmesh);
    runMumps(cmeshMumps, neldiv, pord, nthreads);
  }
  #endif

  return 0;
}
