// main_mumps_3d_complex.cpp
// Demonstrates complex sparse matrix solve using zmumps_c via TPZMumpsSolver<CSTATE>.
// Runs both ZMUMPS and Pardiso (MKL) and compares their solutions.
// Uses TPZL2Projection<CSTATE> to build a Hermitian positive-definite mass matrix.

#include "OutputCapture.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSSpStructMatrixMumps.h"
#include "TPZSYSMPPardiso.h"
#include "TPZSimpleTimer.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzfmatrix.h"
#include "pzlog.h"
#include "pzstepsolver.h"
#include <Projection/TPZL2Projection.h>
#include <TPZGenGrid3D.h>
#include <TPZLinearAnalysis.h>
#include <TPZSYSMPMumps.h>
#include <cmath>
#include <complex>
#include <iostream>

// Forcing function: f(x) = 1 + i (constant complex source)
void ForcingFunctionComplex(const TPZVec<REAL> &pt, TPZVec<CSTATE> &f) {
  f.resize(1);
  f[0] = CSTATE(1.0, 1.0);
}

enum EnumMatids {
  EMatId = 1,
  EBCDirichlet = -1,
  EBCNeumann = -2,
};

TPZGeoMesh *createMesh3D(const TPZVec<int> &nelDiv) {
  TPZVec<REAL> minX(3, 0.);
  TPZVec<REAL> maxX(3, 1.);
  MMeshType elType = MMeshType::ETetrahedral;
  TPZGenGrid3D gen3d(minX, maxX, nelDiv, elType);
  gen3d.BuildVolumetricElements(EMatId);
  // All boundary faces get EBCDirichlet except one face with EBCNeumann
  TPZGeoMesh *gmesh = gen3d.BuildBoundaryElements(EBCDirichlet, EBCNeumann, EBCDirichlet, EBCDirichlet, EBCDirichlet, EBCDirichlet);
  return gmesh;
}

TPZCompMesh *createCompMeshComplex(TPZGeoMesh *gmesh, const int pord) {
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh, /*isComplex=*/true);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetDefaultOrder(pord);
  cmesh->SetAllCreateFunctionsContinuous();

  // TPZL2Projection<CSTATE>: assembles the H1 mass matrix with a complex RHS.
  // This yields a real symmetric positive-definite stiffness matrix (mass matrix)
  // stored as CSTATE, which ZMUMPS can solve.
  auto *mat = new TPZL2Projection<CSTATE>(EMatId, gmesh->Dimension(), 1);
  mat->SetForcingFunction(ForcingFunctionComplex, 2);
  cmesh->InsertMaterialObject(mat);

  // Dirichlet BC: u = 0
  TPZFMatrix<CSTATE> val1(1, 1, 0.);
  TPZManVector<CSTATE, 1> val2(1, CSTATE(0., 0.));
  auto *bcd = mat->CreateBC(mat, EBCDirichlet, 0 /*Dirichlet*/, val1, val2);
  cmesh->InsertMaterialObject(bcd);

  // Neumann BC: flux = 0 (natural, i.e. do nothing — but insert to be explicit)
  auto *bcn = mat->CreateBC(mat, EBCNeumann, 1 /*Neumann*/, val1, val2);
  cmesh->InsertMaterialObject(bcn);

  cmesh->AutoBuild();
  return cmesh;
}

void runPardisoComplex(TPZAutoPointer<TPZCompMesh> cmesh, int neldiv, int pOrder, int nthreads) {
  TPZLinearAnalysis an(cmesh, RenumType::ENone);

  TPZSSpStructMatrix<CSTATE> matsp(cmesh);
  matsp.SetNumThreads(nthreads);
  an.SetStructuralMatrix(matsp);

  TPZStepSolver<CSTATE> step;
  step.SetDirect(ECholesky); // Mass matrix is Hermitian positive definite

  an.SetSolver(step);
  an.Assemble();

  auto *smp = dynamic_cast<TPZSYsmpMatrix<CSTATE> *>(an.MatrixSolver<CSTATE>().Matrix().operator->());
  if (!smp) {
    DebugStop();
    return;
  }
  const int64_t nrows = smp->Rows();
  const int64_t nnz = smp->IA()[nrows];

  auto &mSolverCast = an.MatrixSolver<CSTATE>();
  auto mCast = mSolverCast.Matrix();
  mCast->SetDefPositive(true);
  auto &solverControl = dynamic_cast<TPZSYsmpMatrixPardiso<CSTATE> *>(mCast.operator->())->GetPardisoControl();
  solverControl.SetMessageLevel(1);

  std::cout << "Number of non-zeros: " << nnz << "\n";

  SolverMetrics metrics;
  std::string capturedOutput;
  double neopzSolveTime = 0.0;
  {
    OutputCapture capture;

    TPZSimpleTimer tPZSolve;
    an.Solve();
    neopzSolveTime = tPZSolve.ReturnTimeDouble() / 1000.0;

    capturedOutput = capture.getOutput();
    metrics = capture.parsePardisoOutput(capturedOutput);
  }

  std::cout << "Pardiso Complex Solver Metrics:" << std::endl
            << "nelDiv: " << neldiv << std::endl
            << "pOrder: " << pOrder << std::endl
            << "Number of equations: " << cmesh->NEquations() << std::endl
            << "Estimated NNZ: " << nnz << std::endl;
  if (metrics.numThreads)
    std::cout << "nThreads: " << metrics.numThreads << std::endl;
  if (!metrics.orderingBasedOn.empty())
    std::cout << "Ordering: " << metrics.orderingBasedOn << std::endl;
  if (metrics.factorizationTime)
    std::cout << "factorizationTime: " << metrics.factorizationTime << " s" << std::endl;
  if (metrics.solveTime)
    std::cout << "solveTime: " << std::fixed << std::setprecision(7)
              << metrics.solveTime << " s" << std::endl;
  std::cout << "neopzSolveTime: " << std::fixed << std::setprecision(7)
            << neopzSolveTime << " s" << std::endl;

  an.Solution().Print("Solution (complex) from Pardiso");

  // Post-processing
  const std::string plotfile("PostProcessPardisoComplex");
  constexpr int vtkRes(0);
  TPZManVector<std::string, 1> fields = {"Solution"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.Do();
}

/// Returns ||a - b||_2 / ||a||_2  (or ||a-b||_2 when ||a||_2 ~= 0).
static REAL RelativeL2ErrorComplex(const TPZFMatrix<CSTATE> &ref, const TPZFMatrix<CSTATE> &other) {
  const int64_t n = ref.Rows();
  REAL normDiff = 0., normRef = 0.;
  for (int64_t i = 0; i < n; i++) {
    const CSTATE d = ref.GetVal(i, 0) - other.GetVal(i, 0);
    normDiff += std::norm(d);
    normRef += std::norm(ref.GetVal(i, 0));
  }
  normDiff = std::sqrt(normDiff);
  normRef = std::sqrt(normRef);
  return (normRef > 1e-10) ? normDiff / normRef : normDiff;
}

void runMumpsComplex(TPZAutoPointer<TPZCompMesh> cmesh, int neldiv, int pOrder, int nthreads) {
  TPZLinearAnalysis an(cmesh, RenumType::ENone);

  TPZSSpStructMatrixMumps<CSTATE> matsp(cmesh);
  matsp.SetNumThreads(nthreads);
  an.SetStructuralMatrix(matsp);

  TPZStepSolver<CSTATE> step;
  // Use ELDLt (symmetric indefinite): ZMUMPS has no Hermitian variant.
  // The mass matrix is real-valued stored as CSTATE, so it is both Hermitian
  // and complex-symmetric; SYM=2 (symmetric general) is the correct ZMUMPS mode.
  step.SetDirect(ELDLt);

  an.SetSolver(step);
  an.Assemble();

  auto *smp = dynamic_cast<TPZSYsmpMatrix<CSTATE> *>(an.MatrixSolver<CSTATE>().Matrix().operator->());
  if (!smp) {
    DebugStop();
    return;
  }
  const int64_t nrows = smp->Rows();
  const int64_t nnz = smp->IA()[nrows];

  auto &mSolverCast = an.MatrixSolver<CSTATE>();
  auto mCast = mSolverCast.Matrix();
  auto &solverControl = dynamic_cast<TPZSYsmpMatrixMumps<CSTATE> *>(mCast.operator->())->GetMumpsControl();
  // Explicitly override to SymProp::Sym so that MatrixType() uses SYM=2 (complex symmetric).
  solverControl.SetMatrixType(SymProp::Sym, TPZMumpsSolver<CSTATE>::MProperty::EIndefinite);
  // Lock settings so TPZSYsmpMatrixMumps::Decompose won't re-detect from
  // the matrix's fSymProp (which is Herm for all complex symmetric matrices in NeoPZ).
  solverControl.LockSettings();
  solverControl.SetMessageLevel(2);

  std::cout << "Number of non-zeros: " << nnz << "\n";

  SolverMetrics metrics;
  std::string capturedOutput;
  double neopzSolveTime = 0.0;
  {
    OutputCapture capture;

    TPZSimpleTimer tPZSolve;
    an.Solve();
    neopzSolveTime = tPZSolve.ReturnTimeDouble() / 1000.0;

    capturedOutput = capture.getOutput();
    metrics = capture.parseMumpsOutput(capturedOutput);
  }

  std::cout << "ZMUMPS Complex Solver Metrics:" << std::endl
            << "nelDiv: " << neldiv << std::endl
            << "pOrder: " << pOrder << std::endl
            << "Number of equations: " << cmesh->NEquations() << std::endl
            << "Estimated NNZ: " << nnz << std::endl;
  if (metrics.numThreads)
    std::cout << "nThreads: " << metrics.numThreads << std::endl;
  if (!metrics.orderingBasedOn.empty())
    std::cout << "Ordering: " << metrics.orderingBasedOn << std::endl;
  if (metrics.factorizationTime)
    std::cout << "factorizationTime: " << metrics.factorizationTime << " s" << std::endl;
  if (metrics.solveTime)
    std::cout << "solveTime: " << std::fixed << std::setprecision(7)
              << metrics.solveTime << " s" << std::endl;
  std::cout << "neopzSolveTime: " << std::fixed << std::setprecision(7)
            << neopzSolveTime << " s" << std::endl;

  an.Solution().Print("Solution (complex) from ZMUMPS");

  // Post-processing
  const std::string plotfile("PostProcessMumpsComplex");
  constexpr int vtkRes(0);
  TPZManVector<std::string, 1> fields = {"Solution"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.Do();
}

int main(int argc, char *const argv[]) {
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  const int nthreads = 4;
  int neldiv = 0, n;

  if (argc > 1) {
    neldiv = atoi(argv[1]);
  } else {
    std::cout << "Enter number of divisions: ";
    std::cin >> n;
    neldiv = n;
  }
  std::cout << "Divisions per direction: " << neldiv << std::endl;
  if (neldiv <= 0) {
    std::cerr << "Invalid input, using 2 divisions." << std::endl;
    neldiv = 2;
  }

  TPZGeoMesh *gmesh = createMesh3D({neldiv, neldiv, neldiv});

  const int pOrder = argc > 2 ? atoi(argv[2]) : 1;
  TPZAutoPointer<TPZCompMesh> cmesh = createCompMeshComplex(gmesh, pOrder);

  std::cout << "Number of equations: " << cmesh->NEquations() << "\n";

  TPZFMatrix<CSTATE> solPardiso, solMumps;

#ifdef PZ_USING_MKL
  {
    TPZAutoPointer<TPZCompMesh> cmeshPardiso = new TPZCompMesh(*cmesh);
    runPardisoComplex(cmeshPardiso, neldiv, pOrder, nthreads);
    solPardiso = cmeshPardiso->Solution();
  }
#endif

#ifdef PZ_USING_MUMPS
  {
    TPZAutoPointer<TPZCompMesh> cmeshMumps = new TPZCompMesh(*cmesh);
    runMumpsComplex(cmeshMumps, neldiv, pOrder, nthreads);
    solMumps = cmeshMumps->Solution();
  }
#endif

#ifdef PZ_USING_MKL &&PZ_USING_MUMPS
  const REAL err = RelativeL2ErrorComplex(solPardiso, solMumps);
  std::cout << "\n=== Solution comparison ===" << std::endl
            << "Relative L2 error (Pardiso vs ZMUMPS): " << std::scientific << err << std::endl;
  if (err < 1e-10)
    std::cout << "PASS: solutions agree to within 1e-10." << std::endl;
  else
    std::cout << "WARN: solutions differ more than 1e-10!" << std::endl;
#else
  std::cout << "\n=== Solution comparison skipped (only one solver [or none] available) ===" << std::endl;
#endif

  return 0;
}
