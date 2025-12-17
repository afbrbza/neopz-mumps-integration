#include "TPZMultiphysicsCompMesh.h"
#include "TPZRefPatternDataBase.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSSpStructMatrixMumps.h"
#include "TPZSimpleTimer.h"
#include "TPZSpStructMatrixMumps.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <DarcyFlow/TPZDarcyFlow.h>
#include <TPZGenGrid2D.h>
#include <TPZGenGrid3D.h>
#include <TPZLinearAnalysis.h>
#include <iostream>
#include <pzlog.h>

using namespace std;

enum EnumMatids {
  EMatId = 1, // Material ID for the domain
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
  gmesh->NodeVec()[4] = nod4;
  gmesh->NodeVec()[5] = nod5;

  // ----- First quadrilateral element -----
  TPZManVector<int64_t, 4> nodeIndexes(4);
  nodeIndexes[0] = 0;
  nodeIndexes[1] = 1;
  nodeIndexes[2] = 2;
  nodeIndexes[3] = 3;

  int64_t index = 0;
  TPZGeoEl *gel = gmesh->CreateGeoElement(EQuadrilateral, nodeIndexes, EMatId, index);

  // ----- Second quadrilateral element -----
  nodeIndexes[0] = 1;
  nodeIndexes[1] = 4;
  nodeIndexes[2] = 5;
  nodeIndexes[3] = 2;
  index = 1;
  TPZGeoEl *gel2 = gmesh->CreateGeoElement(EQuadrilateral, nodeIndexes, EMatId, index);

  // ----- Creating BCs -----
  int iddiri = 1;
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

TPZCompMesh *createCompMesh(TPZGeoMesh *gmesh, const int pord) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(pord);
  cmesh->SetAllCreateFunctionsContinuous();

  // Add material (weak formulation)
  TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
  mat->SetConstantPermeability(1.0); // Set constant permeability
  cmesh->InsertMaterialObject(mat);  // insert the material into the computational mesh

  // Add boundary conditions
  int dirichletType = 0 /* Dirichlet condition */, neumannType = 1 /* Neumann condition */;
  int mixedType = 2;                 // Mixed condition
  TPZManVector<REAL, 1> val2(1, 3.); // Parts that goes to the RHS
  TPZFMatrix<REAL> val1(1, 1, 0.);   // Parts that goes to the Stiffness matrix

  // Create and insert boundary conditions
  val2[0] = 0.; // Neumann condition value
  TPZBndCondT<REAL> *bcond =
      mat->CreateBC(mat, EBottom, neumannType, val1, val2); // create BC
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

/**
 * @brief Main function for solving a finite element problem using NeoPZ library.
 *
 * This program demonstrates the complete workflow of solving a PDE using the NeoPZ
 * finite element library:
 * 1. Creates a geometric mesh (either using GenGrid or manual creation)
 * 2. Creates a computational mesh with specified polynomial order
 * 3. Sets up a linear analysis with a chosen direct solver
 * 4. Assembles and solves the linear system
 * 5. Post-processes results to VTK format
 *
 * @details
 * The solver selection includes three options:
 * - EPardiso: Uses Pardiso solver (requires Pardiso installation). The direct solver
 *   type is automatically detected based on the matrix type and the SetDirect() call
 *   is for clarity only to avoid confusion.
 * - EMumps: Uses MUMPS solver (requires MUMPS installation). The direct solver type
 *   is automatically detected based on the matrix type and the SetDirect() call is
 *   for clarity only to avoid confusion.
 * - ESkyline: Uses Skyline solver. The direct solver is explicitly used via SetDirect().
 *
 * @note The Cholesky (LLT) factorization is used as the direct solver. This requires
 * that the matrix is symmetric. For Darcy flow, the matrix is symmetric but INDEFINITE
 * (not positive definite), so MUMPS/Pardiso will use LDL^T factorization with pivoting.
 * Do NOT call SetDefPositive(true) for Darcy flow matrices as they are indefinite.
 *
 * @param argc Number of command-line arguments
 * @param argv Command-line arguments
 * @return int Exit status (0 on success)
 */
int main(int argc, char const *argv[]) {

#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif
  const int nthreads = 0;

  bool print = false;

  // ----- Create geometric mesh -----
  const bool isUseGenGrid = true; // set to 'false' to use manual gmesh creation
  TPZAutoPointer<TPZGeoMesh> gmesh = nullptr;
  if (isUseGenGrid) {
    const int neldiv = 150; // number of elements in each direction
    gmesh = createMeshWithGenGrid({neldiv, neldiv}, {0., 0.}, {2., 1.});
  } else {
    gmesh = createRectangularGmesh();
  }
  if (print) gmesh->Print(std::cout);
  std::ofstream out("geomesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  // ----- Create computational mesh -----
  const int pord = 3; // polynomial order (number of functions per element [number of equations per element]) increases the density of the global matrix
  TPZAutoPointer<TPZCompMesh> cmesh = createCompMesh(gmesh.operator->(), pord);
  if (print) cmesh->Print(std::cout);

  // ----- Create analysis object -----
  TPZLinearAnalysis an(cmesh);

  // criar um enum dos solvers que estou usando: Pardiso, Mumps e o que o Skyline está chamando
  enum EnumSolvers {
    EPardiso,
    EMumps,
    ESkyline
  };

  const EnumSolvers solverType = EMumps;
  TPZAutoPointer<TPZStructMatrixT<STATE>> matsp;
  if (solverType == EPardiso) {
    std::cout << "Using Pardiso solver" << std::endl;
    matsp = new TPZSSpStructMatrix<STATE>(cmesh); /* para quem tem pardiso */
  } else if (solverType == EMumps) {
    // NOTE: Using non-symmetric matrix format for MUMPS due to bugs in symmetric version
    // See BUG_REPORT_MUMPS_SYMMETRIC.md for details
    // std::cout << "Using MUMPS solver (symmetric format)" << std::endl;
    // matsp = new TPZSSpStructMatrixMumps<STATE>(cmesh);
    std::cout << "Using MUMPS solver (non-symmetric format)" << std::endl;
    matsp = new TPZSpStructMatrixMumps<STATE>(cmesh);
  } else {
    std::cout << "Using Skyline solver" << std::endl;
    matsp = new TPZSkylineStructMatrix<STATE>(cmesh);
  }

  matsp->SetNumThreads(nthreads); // number of threads
  an.SetStructuralMatrix(*matsp);
  int64_t neq = cmesh->NEquations();
  std::cout << "Number of equation = " << neq << std::endl;

  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky); // direct solver
  an.SetSolver(step);

  if (solverType == ESkyline) {
    an.Run(); // assembles and solves the linear system
  } else {
    if (neq > 20000)
      TPZSimpleTimer t("Time for assembly", true);
    an.Assemble();

    if (neq > 20000) {
      std::cout << "Entering Solve\n";
      std::cout.flush();
    }

    auto mCast = an.MatrixSolver<STATE>().Matrix();
    mCast->SetDefPositive(true);

    if (neq > 20000)
      TPZSimpleTimer t("Time for solving", true);

    an.Solve();
  }

  print = false;
  if (print) an.Solution().Print("Solution");

  // Post-processing the results
  const std::string plotfile = "postproc"; // without extension
  constexpr int vtkRes(0);                 // 0 for no refinement

  TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

  vtk.Do();

  return 0;
}
