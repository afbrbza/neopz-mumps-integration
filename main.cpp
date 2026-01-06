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
#include <TPZSloanRenumbering.h>
#include <iostream>
#include <mkl_service.h>
#include <numeric>
#include <omp.h>
#include <pzlog.h>
#include <unistd.h>

extern "C" void openblas_set_num_threads(int num_threads);

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

  // Apply Sloan renumbering to reduce bandwidth
  {
    TPZSimpleTimer timer("Sloan Renumbering");
    std::cout << "Applying Sloan renumbering to reduce bandwidth..." << std::endl;

    TPZVec<int64_t> perm, iperm;
    TPZSloanRenumbering sloan;
    sloan.SetElementsNodes(cmesh->NElements(), cmesh->NConnects());

    // Build element graph
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    cmesh->ComputeElGraph(elgraph, elgraphindex);
    sloan.SetElementGraph(elgraph, elgraphindex);

    // Perform reordering
    sloan.Resequence(perm, iperm);

    // Apply permutation to mesh
    cmesh->Permute(perm);

    std::cout << "Renumbering completed in " << timer.ReturnTimeDouble() / 1000.0 << " seconds" << std::endl;
  }

  return cmesh;
}

std::string executeCommand(const std::string &command) {
  std::array<char, 128> buffer;
  std::string result;
  FILE *pipe = popen(command.c_str(), "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
    result += buffer.data();
  }
  int status = pclose(pipe);
  if (status == -1) throw std::runtime_error("pclose() failed!");
  return result;
}

int main(int argc, char const *argv[]) {
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  enum Solvers { EPardiso,
                 EMumps };
  const int nthreadsAssemble = 0; // number of threads (0 = automatic)
  auto auxExecSolvers = [&](Solvers solver,
                            TPZLinearAnalysis &an,
                            double &referenceSolveTime,
                            int nthreads,
                            bool doAssemble = true,
                            double assemblyTime = 0.0)
      -> map<string, double> {
    // var to save result of solver execution an return it
    map<string, double> result;

    if (doAssemble) {
      double assemblyTimeLocal = 0;
      {
        TPZSimpleTimer t("Time for assembly", true);
        an.Assemble();
        assemblyTimeLocal = t.ReturnTimeDouble() / 1000;
        result["AssemblyTime_s"] = assemblyTimeLocal;
      }
      cout << "Seconds for assembly: " << result["AssemblyTime_s"] << " s\n";
    } else {
      result["AssemblyTime_s"] = assemblyTime;
      cout << "Reusing assembly from first iteration\n";
    }

    if (solver == EPardiso) {
      result["RealUsedThreadsAssembly"] = mkl_get_max_threads();
      cout << "MKL threads used in assembly: " << result["RealUsedThreadsAssembly"] << "\n";
    } else if (solver == EMumps) {
#pragma omp parallel
      {
#pragma omp master
        {
          result["RealUsedThreadsAssembly"] = omp_get_num_threads();
          cout << "OpenMP threads used in assembly: " << result["RealUsedThreadsAssembly"] << "\n";
        }
      }
    }

    auto mCast = an.MatrixSolver<STATE>().Matrix();
    mCast->SetDefPositive(true);

    // Solve with timing
    double solveTime = 0;
    {
      TPZSimpleTimer t("Time for solve", true);
      an.Solve();
      solveTime = t.ReturnTimeDouble() / 1000;
      result["SolveTime_s"] = solveTime;
    }
    cout << "Seconds for solve: " << result["SolveTime_s"] << " s\n";

    if (nthreads == 1) {
      referenceSolveTime = result["SolveTime_s"];
      result["SpeedUP"] = 1.0;
    } else {
      result["SpeedUP"] = result["SolveTime_s"] > 0 ? ((referenceSolveTime ?: 0) / result["SolveTime_s"]) : 0;
    }
    cout << "SpeedUP 1/" << nthreads << ": " << result["SpeedUP"] << "\n";

    if (solver == EPardiso) {
      result["RealUsedThreadsSolve"] = mkl_get_max_threads();
      cout << "MKL threads used in solve: " << result["RealUsedThreadsSolve"] << "\n";
    } else if (solver == EMumps) {
#pragma omp parallel
      {
#pragma omp master
        {
          result["RealUsedThreadsSolve"] = omp_get_num_threads();
          cout << "OpenMP threads used in solve: " << result["RealUsedThreadsSolve"] << "\n";
        }
      }
    }

    return result;
  };

  const TPZVec<int> nElemsDiv{250};
  const TPZVec<int> POrds{2};
  /* ************************************************************ */
  // int maxNumOfThreads = executeCommand("nproc --all").empty() ? 8 : std::stoi(executeCommand("nproc --all"));
  // TPZVec<int> nThreadsSolver(maxNumOfThreads);
  // iota(nThreadsSolver.begin(), nThreadsSolver.end(), 1);
  /* ************************************************************ */
  TPZVec<int> nThreadsSolver{1, 12};
  /* ************************************************************ */

  const bool isOutFile = true;
  std::ofstream csvFile("solver_performance_results.csv");
  std::ostream &resultsFile = isOutFile ? (std::ostream &)csvFile : std::cout;
  resultsFile << "NumElementsPerDir,PolyOrder,NumEquations,Solver,NumSolverThreads,AssemblyTime_s,RealUsedThreadsAssembly,SolveTime_s,RealUsedThreadsSolve,SpeedUP\n";

  // executeCommand("clear"); // clear console output
  const bool execPardiso = true;
  const bool execMumps = true;

  double referenceSolveTimePardiso = 0.0;
  double referenceSolveTimeMumps = 0.0;

  for (auto neldiv : nElemsDiv) {
    TPZAutoPointer<TPZGeoMesh> gmesh = createMeshWithGenGrid({neldiv, neldiv}, {0., 0.}, {2., 1.});

    for (auto pord : POrds) {
      // Create the original cmesh once per (neldiv, pord) combination
      TPZAutoPointer<TPZCompMesh> cmesh_original = createCompMesh(gmesh.operator->(), pord);

      int64_t neq = cmesh_original->NEquations();

      map<string, double> result;
      double pardisoAssemblyTime = 0.0;
      if (execPardiso) {
        bool firstIteration = true;

        // Clone cmesh ONCE for all Pardiso iterations
        TPZAutoPointer<TPZCompMesh> cmesh(cmesh_original->Clone());

        // Create analysis object ONCE for all Pardiso iterations
        TPZLinearAnalysis anPardiso(cmesh);
        TPZAutoPointer<TPZStructMatrixT<STATE>> matspPardiso =
            new TPZSSpStructMatrix<STATE>(cmesh);
        matspPardiso->SetNumThreads(nthreadsAssemble);
        anPardiso.SetStructuralMatrix(*matspPardiso);
        TPZStepSolver<STATE> stepPardiso;
        stepPardiso.SetDirect(ECholesky);
        anPardiso.SetSolver(stepPardiso);

        for (auto nthreads : nThreadsSolver) {
          // ===== PARDISO TEST =====
          {
            std::cout << "\n\nRunning for " << neldiv << " elements in each direction.\n";
            std::cout << "  Polynomial order: " << pord << std::endl;
            std::cout << "Number of equation = " << neq << std::endl;
            std::cout << "\n=== Using Pardiso solver with " << nthreads << " threads ===\n";

            // Pardiso uses MKL, disable OpenBLAS to avoid conflicts
            openblas_set_num_threads(1);

            // Disable MKL dynamic threading to force exact number of threads
            mkl_set_dynamic(0);

            // Set number of threads for MKL and OpenMP
            omp_set_num_threads(nthreads);
            mkl_set_num_threads(nthreads);

            // Also set domain-specific threads
            mkl_domain_set_num_threads(nthreads, MKL_DOMAIN_ALL);

            // Verify MKL settings
            int mkl_max_threads = mkl_get_max_threads();
            std::cout << "MKL max threads: " << mkl_max_threads << std::endl;

            result = auxExecSolvers(EPardiso, anPardiso, referenceSolveTimePardiso, nthreads, firstIteration, pardisoAssemblyTime);
            if (firstIteration) {
              pardisoAssemblyTime = result["AssemblyTime_s"];
              firstIteration = false;
            }
          }
          resultsFile << neldiv << "," << pord << "," << neq << ",Pardiso," << nthreads << ","
                      << result["AssemblyTime_s"] << "," << result["RealUsedThreadsAssembly"] << ","
                      << result["SolveTime_s"] << "," << result["RealUsedThreadsSolve"]
                      << "," << result["SpeedUP"] << "\n";
        }
      }

      if (execMumps) {
        bool firstIteration = true;
        double mumpsAssemblyTime = 0.0;

        // Clone cmesh ONCE for all MUMPS iterations
        TPZAutoPointer<TPZCompMesh> cmesh(cmesh_original->Clone());

        // Create analysis object ONCE for all MUMPS iterations
        TPZLinearAnalysis anMumps(cmesh);
        TPZAutoPointer<TPZStructMatrixT<STATE>> matspMumps =
            new TPZSpStructMatrixMumps<STATE>(cmesh);
        matspMumps->SetNumThreads(nthreadsAssemble);
        anMumps.SetStructuralMatrix(*matspMumps);
        TPZStepSolver<STATE> stepMumps;
        stepMumps.SetDirect(ECholesky);
        anMumps.SetSolver(stepMumps);

        for (auto nthreads : nThreadsSolver) {
          // ===== MUMPS TEST =====
          {
            std::cout << "\n\nRunning for " << neldiv << " elements in each direction.\n";
            std::cout << "  Polynomial order: " << pord << std::endl;
            std::cout << "Number of equation = " << neq << std::endl;
            std::cout << "\n=== Using MUMPS solver with " << nthreads << " threads ===\n";

            // MUMPS uses OpenBLAS for BLAS operations + OpenMP for parallelism
            mkl_set_num_threads(1); // Disable MKL to avoid conflicts
            openblas_set_num_threads(nthreads);
            omp_set_num_threads(nthreads);

// Verify OpenMP settings
#pragma omp parallel
            {
#pragma omp master
              {
                std::cout << "OpenMP threads: " << omp_get_num_threads() << std::endl;
              }
            }

            result = auxExecSolvers(EMumps, anMumps, referenceSolveTimeMumps, nthreads, firstIteration, mumpsAssemblyTime);
            if (firstIteration) {
              mumpsAssemblyTime = result["AssemblyTime_s"];
              firstIteration = false;
            }
          }
          resultsFile << neldiv << "," << pord << "," << neq << ",MUMPS," << nthreads << ","
                      << result["AssemblyTime_s"] << "," << result["RealUsedThreadsAssembly"] << ","
                      << result["SolveTime_s"] << "," << result["RealUsedThreadsSolve"]
                      << "," << result["SpeedUP"]
                      << "\n";
        }
      }
    }
  }

  if (isOutFile)
    csvFile.close();
}

// /**
//  * @brief Main function for solving a finite element problem using NeoPZ library.
//  *
//  * This program demonstrates the complete workflow of solving a PDE using the NeoPZ
//  * finite element library:
//  * 1. Creates a geometric mesh (either using GenGrid or manual creation)
//  * 2. Creates a computational mesh with specified polynomial order
//  * 3. Sets up a linear analysis with a chosen direct solver
//  * 4. Assembles and solves the linear system
//  * 5. Post-processes results to VTK format
//  *
//  * @details
//  * The solver selection includes three options:
//  * - EPardiso: Uses Pardiso solver (requires Pardiso installation). The direct solver
//  *   type is automatically detected based on the matrix type and the SetDirect() call
//  *   is for clarity only to avoid confusion.
//  * - EMumps: Uses MUMPS solver (requires MUMPS installation). The direct solver type
//  *   is automatically detected based on the matrix type and the SetDirect() call is
//  *   for clarity only to avoid confusion.
//  * - ESkyline: Uses Skyline solver. The direct solver is explicitly used via SetDirect().
//  *
//  * @note The Cholesky (LLT) factorization is used as the direct solver. This requires
//  * that the matrix is symmetric. For Darcy flow, the matrix is symmetric but INDEFINITE
//  * (not positive definite), so MUMPS/Pardiso will use LDL^T factorization with pivoting.
//  * Do NOT call SetDefPositive(true) for Darcy flow matrices as they are indefinite.
//  *
//  * @param argc Number of command-line arguments
//  * @param argv Command-line arguments
//  * @return int Exit status (0 on success)
//  */
// int mainOriginal(int argc, char const *argv[]) {

// #ifdef PZ_LOG
//   TPZLogger::InitializePZLOG();
// #endif
//   const int nthreads = 0;

//   bool print = false;

//   // ----- Create geometric mesh -----
//   const bool isUseGenGrid = true; // set to 'false' to use manual gmesh creation
//   TPZAutoPointer<TPZGeoMesh> gmesh = nullptr;
//   if (isUseGenGrid) {
//     const int neldiv = 2; // number of elements in each direction
//     gmesh = createMeshWithGenGrid({neldiv, neldiv}, {0., 0.}, {2., 1.});
//   } else {
//     gmesh = createRectangularGmesh();
//   }
//   if (print) gmesh->Print(std::cout);
//   std::ofstream out("geomesh.vtk");
//   TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

//   // ----- Create computational mesh -----
//   const int pord = 1; // polynomial order (number of functions per element [number of equations per element]) increases the density of the global matrix
//   TPZAutoPointer<TPZCompMesh> cmesh = createCompMesh(gmesh.operator->(), pord);
//   if (print) cmesh->Print(std::cout);

//   // ----- Create analysis object -----
//   TPZLinearAnalysis an(cmesh);

//   // criar um enum dos solvers que estou usando: Pardiso, Mumps e o que o Skyline está chamando
//   enum EnumSolvers {
//     EPardiso,
//     EMumps,
//     ESkyline
//   };

//   const EnumSolvers solverType = EMumps;
//   TPZAutoPointer<TPZStructMatrixT<STATE>> matsp;
//   if (solverType == EPardiso) {
//     std::cout << "Using Pardiso solver" << std::endl;
//     matsp = new TPZSSpStructMatrix<STATE>(cmesh); /* para quem tem pardiso */
//   } else if (solverType == EMumps) {
//     // NOTE: Using non-symmetric matrix format for MUMPS due to bugs in symmetric version
//     // See BUG_REPORT_MUMPS_SYMMETRIC.md for details
//     // std::cout << "Using MUMPS solver (symmetric format)" << std::endl;
//     // matsp = new TPZSSpStructMatrixMumps<STATE>(cmesh);
//     std::cout << "Using MUMPS solver (non-symmetric format)" << std::endl;
//     matsp = new TPZSpStructMatrixMumps<STATE>(cmesh);
//   } else {
//     std::cout << "Using Skyline solver" << std::endl;
//     matsp = new TPZSkylineStructMatrix<STATE>(cmesh);
//   }

//   matsp->SetNumThreads(nthreads); // number of threads
//   an.SetStructuralMatrix(*matsp);
//   int64_t neq = cmesh->NEquations();
//   std::cout << "Number of equation = " << neq << std::endl;

//   TPZStepSolver<STATE> step;
//   step.SetDirect(ECholesky); // direct solver
//   an.SetSolver(step);

//   if (solverType == ESkyline) {
//     an.Run(); // assembles and solves the linear system
//   } else {
//     if (neq > 20000) {
//       TPZSimpleTimer t("Time for assembly", true);
//       an.Assemble();
//     } else {
//       an.Assemble();
//     }

//     if (neq > 20000) {
//       std::cout << "Entering Solve\n";
//       std::cout.flush();
//     }

//     auto mCast = an.MatrixSolver<STATE>().Matrix();
//     mCast->SetDefPositive(true);

//     if (neq > 20000) {
//       TPZSimpleTimer t("Time for solve", true);
//       an.Solve();
//     } else {
//       an.Solve();
//     }
//   }

//   print = false;
//   if (print) an.Solution().Print("Solution");

//   // Post-processing the results
//   const std::string plotfile = "postproc"; // without extension
//   constexpr int vtkRes(0);                 // 0 for no refinement

//   TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
//   auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

//   vtk.Do();

//   return 0;
// }
