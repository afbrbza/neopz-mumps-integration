#include "OutputCapture.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrixMumps.h"
#include "TPZSimpleTimer.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzfmatrix.h"
#include "pzlog.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <DarcyFlow/TPZDarcyFlow.h>
#include <TPZGenGrid2D.h>
#include <TPZGenGrid3D.h>
#include <TPZLinearAnalysis.h>
#include <TPZSYSMPMumps.h>
#include <TPZSloanRenumbering.h>
#include <iostream>
#include <omp.h>

using namespace std;

enum EnumMatids {
  EMatId = 1,
  EBottom = 2,
  ETop = 3,
  ELeft = 4,
  ERight = 5
};

TPZGeoMesh *createMeshWithGenGrid(const TPZVec<int> &nelDiv, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX) {
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  TPZGenGrid2D generator(nelDiv, minX, maxX);
  generator.SetElementType(MMeshType::EQuadrilateral);
  generator.Read(gmesh, EMatId);
  generator.SetBC(gmesh, 4, EBottom); // bottom
  generator.SetBC(gmesh, 5, ERight);  // right
  generator.SetBC(gmesh, 6, ETop);    // top
  generator.SetBC(gmesh, 7, ELeft);   // left

  return gmesh;
}

TPZGeoMesh *createRectangularGmesh() {
  TPZGeoMesh *gmesh = new TPZGeoMesh;

  TPZManVector<REAL, 3> coord(3, 0.);

  TPZGeoNode nod0(0, coord, *gmesh); // Node at (0,0,0)
  coord[0] = 1.;
  TPZGeoNode nod1(1, coord, *gmesh); // Node at (1,0,0)
  coord[1] = 1.;
  TPZGeoNode nod2(2, coord, *gmesh); // Node at (1,1,0)
  coord[0] = 0.;
  TPZGeoNode nod3(3, coord, *gmesh); // Node at (0,1,0)
  coord[0] = 2.;
  coord[1] = 0.;
  TPZGeoNode nod4(4, coord, *gmesh); // Node at (2,0,0)
  coord[1] = 1.;
  TPZGeoNode nod5(5, coord, *gmesh); // Node at (2,1,0)

  gmesh->NodeVec().Resize(6);
  gmesh->NodeVec()[0] = nod0;
  gmesh->NodeVec()[1] = nod1;
  gmesh->NodeVec()[2] = nod2;
  gmesh->NodeVec()[3] = nod3;
  gmesh->NodeVec()[3] = nod4;
  gmesh->NodeVec()[3] = nod5;

  // ----- First quadrilateral element -----
  TPZManVector<int64_t, 4> nodeindexes(4);
  nodeindexes[0] = 0;
  nodeindexes[1] = 1;
  nodeindexes[2] = 2;
  nodeindexes[3] = 3;

  int64_t index = 0;
  TPZGeoEl *gel = gmesh->CreateGeoElement(EQuadrilateral, nodeindexes, EMatId, index);

  // ----- Second quadrilateral element -----
  nodeindexes[0] = 1;
  nodeindexes[1] = 4;
  nodeindexes[2] = 5;
  nodeindexes[3] = 2;
  index = 1;
  TPZGeoEl *gel2 = gmesh->CreateGeoElement(EQuadrilateral, nodeindexes, EMatId, index);

  // ----- Creating BCs -----
  gel->CreateBCGeoEl(4, EBottom);
  gel2->CreateBCGeoEl(4, EBottom);
  gel2->CreateBCGeoEl(5, ERight);
  gel2->CreateBCGeoEl(6, ETop);
  gel->CreateBCGeoEl(6, ETop);
  gel->CreateBCGeoEl(7, ELeft);

  gmesh->SetDimension(2);
  gmesh->BuildConnectivity();

  return gmesh;
}

TPZCompMesh *createCompMesh(TPZGeoMesh *gmesh, int pord = 1) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(pord);
  cmesh->SetAllCreateFunctionsContinuous();

  // Add material (weak formulation)
  TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
  mat->SetConstantPermeability(1.0); // Set constant permeability
  cmesh->InsertMaterialObject(mat);  // insert the material into the computational mesh

  // Add boundary conditions
  int dirichletType = 0;             // Dirichlet condition
  int neumannType = 1;               // Neumann condition
  int mixedType = 2;                 // Mixed condition
  TPZManVector<REAL, 1> val2(1, 3.); // Parts that goes to the RHS
  TPZFMatrix<REAL> val1(1, 1, 0.);   // Parts that goes to the Stiffness matrix

  // Create and insert boundary conditions
  val2[0] = 0.;                                                                    // Neumann condition value
  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, neumannType, val1, val2); // create BC
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 1.;                                                  // Dirichlet condition value
  bcond = mat->CreateBC(mat, ERight, dirichletType, val1, val2); // create BC
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 0.;                                              // Dirichlet condition value
  bcond = mat->CreateBC(mat, ETop, neumannType, val1, val2); // create BC
  cmesh->InsertMaterialObject(bcond);

  val2[0] = 3.;                                                 // Dirichlet condition value
  bcond = mat->CreateBC(mat, ELeft, dirichletType, val1, val2); // create BC
  cmesh->InsertMaterialObject(bcond);

  cmesh->AutoBuild(); // create the computational elements and connects

  return cmesh;
}

int main(int argc, char const *argv[]) {
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  TPZSimpleTimer totalTime;

  // nthreads use in NeoPZ assembly
  const int nthreadsAssembly = 0;
  const bool print = false;

  int const num_procs = argc > 1 ? atoi(argv[1]) : 0;
  if (num_procs > 0) {
    omp_set_num_threads(num_procs);
  }

  // ----- Create geometric mesh -----
  const bool isUseGenGrid = true; // set to 'false' to use manual gmesh creation
  const int neldiv = argc > 2 ? atoi(argv[2]) : 140;
  TPZGeoMesh *gmesh = createMeshWithGenGrid({neldiv, neldiv}, {0., 0.}, {2., 1.});

  if (print)
    gmesh->Print(std::cout);

  if (print) {
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  // ----- Create computational mesh -----
  const int pOrder = argc > 3 ? atoi(argv[3]) : 7;
  TPZCompMesh *cmesh = createCompMesh(gmesh, pOrder);
  if (print)
    cmesh->Print(std::cout);

  // ----- Create analysis object -----
  TPZLinearAnalysis an(cmesh);
  // ---- para quem tem mumps
  TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOR<STATE>> matsp(cmesh);
  // TPZSkylineStructMatrix<STATE> matsp(cmesh);
  matsp.SetNumThreads(nthreadsAssembly); // number of threads
  an.SetStructuralMatrix(matsp);
  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky); // direct solver
  an.SetSolver(step);

  auto stiff = an.StructMatrix();
  auto sparseMatrix = dynamic_cast<TPZSYsmpMatrix<STATE> *>(stiff->Create());
  MUMPS_INT nnz = sparseMatrix->A().size();

  an.Assemble();

  // auto & means that mSolverCast is a reference to the object returned by MatrixSolver<STATE>()
  auto &mSolverCast = an.MatrixSolver<STATE>();

  // auto without & means that mCast is a copy of the smart pointer returned by Matrix()
  auto mCast = mSolverCast.Matrix();
  mCast->SetDefPositive(true);

  // auto & here means that mumpsControl is a reference to the object returned by GetMumpsControl()
  // mCast.operator->() returns a pointer to the matrix being used in the solver
  auto &mumpsControl = dynamic_cast<TPZSYsmpMatrixMumps<STATE> *>(mCast.operator->())->GetMumpsControl();

  mumpsControl.SetMessageLevel(2); // ICNTL(4) = 2 (detailed statistics)

  SolverMetrics metrics;
  std::string capturedOutput;
  double neopzSolveTime = 0.0;
  {
    OutputCapture capture;

    TPZSimpleTimer t;
    an.Solve();
    neopzSolveTime = t.ReturnTimeDouble() / 1000.0; // convert to seconds

    capturedOutput = capture.getOutput();
    metrics = capture.parseMumpsOutput(capturedOutput);
  }

  // std::cout << capturedOutput << std::endl;

  std::cout << "MUMPS Solver Metrics:" << std::endl
            << "nelDiv: " << neldiv << std::endl
            << "pOrder: " << pOrder << std::endl
            << "Number of equations: " << cmesh->NEquations() << std::endl
            << "Estimated number of non-zeros: " << nnz << std::endl
            << "nThreads: " << metrics.numThreads << std::endl
            << "factorizationTime: " << metrics.factorizationTime << " s" << std::endl
            << "solveTime: " << std::fixed << std::setprecision(7) << metrics.solveTime << " s" << std::endl
            << "neopzSolveTime: " << std::fixed << std::setprecision(7) << neopzSolveTime << " s" << std::endl;

  // an.Run(); // assembles and solves the linear system
    an.Solution().Print("Solution");

  if (print) {
    // Post-processing the results
    const std::string plotfile = "postproc"; // without extension
    constexpr int vtkRes(0);                 // 0 for no refinement

    TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

    vtk.Do();
  }

  std::cout << std::endl
            << "Total execution time: " << totalTime.ReturnTimeDouble() / 1000.0 << " s" << std::endl;

  return 0;
}
