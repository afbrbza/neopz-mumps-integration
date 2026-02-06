#include "OutputCapture.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
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
#include "pzskylstrmatrix.h"
#include <TPZLinearAnalysis.h>
#include "pzskylmat.h"
#include <TPZSYSMPMumps.h>
#include <TPZSloanRenumbering.h>
#include <iostream>
#include <filesystem>

// Simple constant forcing function for the domain
void ForcingFunction(const TPZVec<REAL> &pt, TPZVec<STATE> &f)
{
  f.resize(1);
  f[0] = 2.0; // Constant source term
}

enum EnumMatids
{
  EMatId = 1,
  EBottom = 2,
  ERight = 3,
  ETop = 4,
  ELeft = 5
};

TPZGeoMesh *createMeshWithGenGrid(const TPZVec<int> &nelDiv, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX)
{ // ----------------O QUE É O GENGRID?----------------
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

TPZGeoMesh *createMesh3D(const TPZVec<int> &nelDiv)
{
  TPZVec<REAL> minX(3, 0.);
  TPZVec<REAL> maxX(3, 1.);
  MMeshType elType = MMeshType::ETetrahedral;
  TPZGenGrid3D gen3d(minX, maxX, nelDiv, elType);
  gen3d.BuildVolumetricElements(1);
  TPZGeoMesh *gmesh = gen3d.BuildBoundaryElements(-1, -3, -1, -2, -1, -1);
  return gmesh;
}

TPZCompMesh *createCompMesh3D(TPZGeoMesh *gmesh, const int pord)
{
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

TPZGeoMesh *createRectangularGmesh()
{
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

TPZCompMesh *createCompMesh(TPZGeoMesh *gmesh, const int pord)
{
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

// --- FUNÇÃO ADICIONADA PARA EXPORTAR MATRIZ NO FORMATO I J VAL ---
void ExportMatrixToCoordinateFormat(TPZMatrix<STATE> *mat, std::string filename)
{
  if (!mat)
    return;

  // 1) Try efficient exports by detecting concrete matrix types
  // 1.a) TPZSYsmpMatrix (symmetric sparse Yale-like storage)
  if (auto *smp = dynamic_cast<TPZSYsmpMatrix<STATE> *>(mat))
  {
    const auto &IA = smp->IA();
    const auto &JA = smp->JA();
    const auto &A = smp->A();
    const int64_t n = smp->Rows();
    const int64_t nnz = IA[n];

    std::ofstream out(filename);
    out << n << " " << n << " " << nnz << "\n";
    out.setf(std::ios::scientific);
    out.precision(16);

    for (int64_t row = 0; row < n; ++row)
    {
      for (int64_t k = IA[row]; k < IA[row + 1]; ++k)
      {
        int64_t col = JA[k];
        STATE val = A[k];
        out << row << " " << col << " " << val << "\n";
        // If you want full symmetric matrix (both triangles), optionally emit (col,row,val) when col!=row
      }
    }
    out.close();
    std::cout << "Exportado (TPZSYsmpMatrix) para " << filename << " nnz=" << nnz << "\n";
    double esparcidade = static_cast<double>(100. * nnz / (n * n));
    std::cout << "Esparcidade = " << esparcidade << " % \n";
    return;
  }

  // 1.c) TPZSkylMatrix (skyline storage) - produce output row-by-row
  if (auto *sky = dynamic_cast<TPZSkylMatrix<STATE> *>(mat))
  {
    const int64_t n = sky->Rows();
    // First compute nnz by summing skyline heights + diagonal
    int64_t nnz = 0;
    for (int64_t col = 0; col < n; ++col)
    {
      int64_t height = sky->SkyHeight(col);
      nnz += (height + 1); // inclusive of diagonal
    }

    std::ofstream out(filename);
    out << n << " " << n << " " << nnz << "\n";
    out.setf(std::ios::scientific);
    out.precision(16);

    // Collect entries grouped by row so we can write all entries of row 0, then row 1, ...
    std::vector<std::vector<std::pair<int64_t, STATE>>> rows;
    rows.resize(n);
    for (int64_t col = 0; col < n; ++col)
    {
      int64_t height = sky->SkyHeight(col);
      int64_t row0 = col - height;
      for (int64_t row = row0; row <= col; ++row)
      {
        STATE val = sky->GetVal(row, col);
        if (std::abs(val) > 1e-16)
          rows[row].emplace_back(col, val);
      }
    }

    for (int64_t row = 0; row < n; ++row)
    {
      for (auto &pr : rows[row])
      {
        out << row << " " << pr.first << " " << pr.second << "\n";
      }
    }

    out.close();
    std::cout << "Exportado (TPZSkylMatrix) para " << filename << " nnz=" << nnz << "\n";
    double esparcidade = static_cast<double>(100. * nnz / (n * n));
    std::cout << "Esparcidade = " << esparcidade << " % \n";
    return;
  }
}
// -----------------------------------------------------------------

// --- FUNÇÃO PARA EXPORTAR VETOR (LOAD VECTOR) ---
void ExportVectorToFile(const TPZFMatrix<STATE> &vec, std::string filename)
{
  int64_t n = vec.Rows();
  std::ofstream out(filename);
  out << n << " " << 1 << "\n";
  out.setf(std::ios::scientific);
  out.precision(16);
  for (int64_t i = 0; i < n; ++i)
  {
    out << vec(i, 0) << "\n";
  }
  out.close();
  std::cout << "Vetor exportado para " << filename << " com " << n << " entradas.\n";
}
// -----------------------------------------------------------------

int main(int argc, char *const argv[])
{
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  static const std::vector<std::string> COMPILE_TIME_OUTPUT_DIR = {
      std::filesystem::path(__FILE__).parent_path().parent_path().string() + "/build/output/"};

  std::vector<std::filesystem::path> outPaths;
  if (!COMPILE_TIME_OUTPUT_DIR.empty())
  {
    for (const auto &d : COMPILE_TIME_OUTPUT_DIR)
    {
      std::filesystem::path p(d);
      std::error_code ec;
      std::filesystem::create_directories(p, ec);
      if (ec)
      {
        std::cerr << "Aviso: nao foi possivel criar o diretorio de saida '" << d << "': " << ec.message() << "\n";
      }
      else
      {
        // std::cout << "Usando diretorio de saida: " << p << "\n";
      }
      outPaths.push_back(p);
    }
  }
  else
  {
    std::cerr << "Aviso: nenhum diretorio de saida compilado. Os arquivos de saida serao escritos no diretorio atual.\n";
    outPaths.push_back(std::filesystem::current_path());
  }

  //---------------------------Create geometric mesh---------------------------
  const int nthreads = 16;
  const bool isUseGenGrid = true;
  TPZGeoMesh *gmesh = nullptr;
  int neldiv = 0, n;
  const bool is3d = true;
  if (isUseGenGrid)
  {
    if (argc > 1)
    {
      neldiv = atoi(argv[1]);
    }
    else
    {
      std::cout << "Enter number of divisions: ";
      std::cin >> n;
      // neldiv = (std::sqrt(n) - 1);
      neldiv = n;
    }
    std::cout << "Number of divisions in each direction: " << neldiv << std::endl;

    // std::cout<< neldiv;
    std::cout << "Total number of nodes: " << ((neldiv + 1) * (neldiv + 1)) << std::endl;
    if (neldiv <= 0)
    {
      std::cerr << "Invalid input, using 1 division." << std::endl;
      neldiv = 1;
    }
    if (is3d)
    {
      // neldiv = 5;
      gmesh = createMesh3D({neldiv, neldiv, neldiv});
    }
    else
      gmesh = createMeshWithGenGrid({neldiv, neldiv}, {0., 0.}, {1., 1.});
  }
  else
  {
    gmesh = createRectangularGmesh();
  }

  // gmesh->Print(std::cout);
  // std::ofstream out("geomesh.vtk");
  // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  //------------------------Create computational mesh-------------------------
  const int pOrder = argc > 2 ? atoi(argv[2]) : 1;
  const int pord = pOrder; // Polynomial order
  TPZAutoPointer<TPZCompMesh> cmesh = nullptr;
  if (is3d)
  {
    cmesh = createCompMesh3D(gmesh, pord);
  }
  else
  {
    cmesh = createCompMesh(gmesh, pord);
  }

  std::cout << "Number of equations: " << cmesh->NEquations() << "\n";

  // TPZAutoPointer<TPZCompMesh> cmesh = createCompMesh(gmesh, pord);
  // cmesh->Print(std::cout);

  //-------------------------Create analysis object--------------------------
  TPZSimpleTimer tot("total time: ");
  TPZLinearAnalysis an(cmesh, RenumType::ENone);
  TPZSSpStructMatrix<STATE> matsp(cmesh);
  // TPZSkylineStructMatrix<STATE> matsp(cmesh);
  matsp.SetNumThreads(nthreads);
  an.SetStructuralMatrix(matsp);

  // Suppress frontal assembly messages
  // std::ofstream devnull("/dev/null");
  // std::streambuf* old_cout = std::cout.rdbuf(devnull.rdbuf());
  // an.Assemble();
  // auto &matrixSolver = dynamic_cast<>MatrixSolver<STATE>().GetPardisoControl();
  // matrixSolver.Mess
  // std::cout.rdbuf(old_cout);

  // TPZAutoPointer<TPZMatrix<STATE>> matrizPreenchida = an.MatrixSolver<STATE>().Matrix();
  // matrizPreenchida->Print("Matriz preenchida:");

  //----------------------PARA RESOLVER A MATRIZ COM PARDISO----------------------
  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky);

  an.SetSolver(step);

  // auto stiff = an.StructMatrix();
  // auto sparseMatrix = dynamic_cast<TPZSYsmpMatrix<STATE> *>(stiff->Create());
  // auto mCast = mSolverCast.Matrix();
  // mCast.A();
  // int64_t nnz = sparseMatrix->A().size();

  an.Assemble();

  // matspauto stiff = an.StructMatrix();
  auto *smp = dynamic_cast<TPZSYsmpMatrix<STATE> *>(an.MatrixSolver<STATE>().Matrix().operator->());
  if (!smp)
  {
    DebugStop();
    return -1;
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

  // auto &mSolverCast = an.MatrixSolver<STATE>();
  // auto mCast = mSolverCast.Matrix();
  // auto &mumpsControl = dynamic_cast<TPZSYsmpMatrixMumps<STATE> *>(mCast.operator->())->GetMumpsControl();
  // mumpsControl.SetMessageLevel(2);

  std::cout << "Number of non-zeros: " << nnz << "\n";
  TPZSimpleTimer t("Solving system");

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
  std::cout << "Tempo de resolução: " << t.ReturnTimeDouble() / 1000. << " segundos.\n";

  // an.Run();
  // an.Solution().Print("Solution");

  // --- CHAMADA DA FUNÇÃO DE EXPORTAÇÃO (em todos os diretórios) ---

  /*for (const auto &outPath : outPaths)
  {
    std::string matrixFilename = (outPath / "LHS.txt").string();
    ExportMatrixToCoordinateFormat(matrizPreenchida.operator->(), matrixFilename);

    std::string loadVectorFilename = (outPath / "RHS.txt").string();
    ExportVectorToFile(an.Rhs(), loadVectorFilename);
  }*/

  // ---------------------------------------

  // Post-processing
  // const std::string plotfile("PostProcess");
  // constexpr int vtkRes(0);

  // TPZManVector<std::string,2> fields = {"Flux","Pressure"};
  // auto vtk = TPZVTKGenerator(cmesh,fields,plotfile, vtkRes);
  // vtk.Do();
  std::cout << "Tempo de resolução total: " << tot.ReturnTimeDouble() / 1000. << " segundos.\n";

  return 0;
}
