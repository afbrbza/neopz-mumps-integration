#include "TPZSpoolesSolver.h"
#include "TPZYSMPSpooles.h"
#include "pzlog.h"
#include "pzmatrix.h"

extern "C" {
#include "misc.h"
#include "InpMtx.h"
#include "ETree.h"
#include "IVL.h"
#include "FrontMtx.h"
#include "SubMtxManager.h"
#include "DenseMtx.h"
#include "SymbFac.h"
#include "Graph.h"
}

#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.spoolescontrol");
#endif

template <class TVar>
TPZSpoolesSolver<TVar>::TPZSpoolesSolver() : TPZMatrixSolver<TVar>() {}

template <class TVar>
TPZSpoolesSolver<TVar>::TPZSpoolesSolver(TPZSpoolesSolver &&copy) noexcept 
    : TPZMatrixSolver<TVar>(std::move(copy)),
      fSymmetry(copy.fSymmetry),
      fProperty(copy.fProperty),
      fInputMatrix(copy.fInputMatrix),
      fETree(copy.fETree),
      fSymbFac(copy.fSymbFac),
      fFrontMtx(copy.fFrontMtx),
      fMtxManager(copy.fMtxManager),
      fSolutionMtx(copy.fSolutionMtx),
      fRhsMtx(copy.fRhsMtx),
      fOldToNewIV(copy.fOldToNewIV),
      fNewToOldIV(copy.fNewToOldIV),
      fMessageLevel(copy.fMessageLevel),
      fNumThreads(copy.fNumThreads),
      fMatrixType(copy.fMatrixType),
      fDecomposed(copy.fDecomposed),
      fSPOOLESInitialized(copy.fSPOOLESInitialized) {
  // Mark source as uninitialized so it won't try to free SPOOLES memory
  copy.fInputMatrix = nullptr;
  copy.fETree = nullptr;
  copy.fSymbFac = nullptr;
  copy.fFrontMtx = nullptr;
  copy.fMtxManager = nullptr;
  copy.fSolutionMtx = nullptr;
  copy.fRhsMtx = nullptr;
  copy.fOldToNewIV = nullptr;
  copy.fNewToOldIV = nullptr;
  copy.fSPOOLESInitialized = false;
  copy.fDecomposed = false;
}

template <class TVar>
TPZSpoolesSolver<TVar>& TPZSpoolesSolver<TVar>::operator=(TPZSpoolesSolver &&copy) noexcept {
  if (this != &copy) {
    FreeSPOOLESMemory();
    
    TPZMatrixSolver<TVar>::operator=(std::move(copy));
    
    fSymmetry = copy.fSymmetry;
    fProperty = copy.fProperty;
    fInputMatrix = copy.fInputMatrix;
    fETree = copy.fETree;
    fSymbFac = copy.fSymbFac;
    fFrontMtx = copy.fFrontMtx;
    fMtxManager = copy.fMtxManager;
    fSolutionMtx = copy.fSolutionMtx;
    fRhsMtx = copy.fRhsMtx;
    fOldToNewIV = copy.fOldToNewIV;
    fNewToOldIV = copy.fNewToOldIV;
    fMessageLevel = copy.fMessageLevel;
    fNumThreads = copy.fNumThreads;
    fMatrixType = copy.fMatrixType;
    fDecomposed = copy.fDecomposed;
    fSPOOLESInitialized = copy.fSPOOLESInitialized;
    
    copy.fInputMatrix = nullptr;
    copy.fETree = nullptr;
    copy.fSymbFac = nullptr;
    copy.fFrontMtx = nullptr;
    copy.fMtxManager = nullptr;
    copy.fSolutionMtx = nullptr;
    copy.fRhsMtx = nullptr;
    copy.fOldToNewIV = nullptr;
    copy.fNewToOldIV = nullptr;
    copy.fSPOOLESInitialized = false;
    copy.fDecomposed = false;
  }
  return *this;
}

template <class TVar>
void TPZSpoolesSolver<TVar>::FreeSPOOLESMemory() {
  if (fSPOOLESInitialized) {
    if (fInputMatrix) InpMtx_free(fInputMatrix);
    if (fETree) ETree_free(fETree);
    if (fSymbFac) IVL_free(fSymbFac);
    if (fFrontMtx) FrontMtx_free(fFrontMtx);
    if (fMtxManager) SubMtxManager_free(fMtxManager);
    if (fSolutionMtx) DenseMtx_free(fSolutionMtx);
    if (fRhsMtx) DenseMtx_free(fRhsMtx);
    if (fOldToNewIV) IV_free(fOldToNewIV);
    if (fNewToOldIV) IV_free(fNewToOldIV);
    
    fInputMatrix = nullptr;
    fETree = nullptr;
    fSymbFac = nullptr;
    fFrontMtx = nullptr;
    fMtxManager = nullptr;
    fSolutionMtx = nullptr;
    fRhsMtx = nullptr;
    fOldToNewIV = nullptr;
    fNewToOldIV = nullptr;
    fSPOOLESInitialized = false;
    fDecomposed = false;
  }
}

template <class TVar>
TPZSpoolesSolver<TVar>::~TPZSpoolesSolver() {
  FreeSPOOLESMemory();
}

template <class TVar>
void TPZSpoolesSolver<TVar>::SetMatrix(TPZAutoPointer<TPZBaseMatrix> Refmat) {
  TPZMatrixSolver<TVar>::SetMatrix(Refmat);
  fDecomposed = false;
}

template <class TVar>
void TPZSpoolesSolver<TVar>::SetMatrixType(SymProp symtype, MProperty prop) {
  fSymmetry = symtype;
  fProperty = prop;
}

template <class TVar>
void TPZSpoolesSolver<TVar>::SetMessageLevel(int lvl) {
  fMessageLevel = lvl;
}

template <class TVar>
TPZSpoolesSolver<TVar> *TPZSpoolesSolver<TVar>::Clone() const {
  std::cerr << "TPZSpoolesSolver::Clone() - Cannot clone SPOOLES solver\n";
  DebugStop();
  return nullptr;
}

template <class TVar>
void TPZSpoolesSolver<TVar>::Decompose() {
  if (fDecomposed) {
    std::cout << "TPZSpoolesSolver::Decompose() - Matrix already decomposed\n";
    return;
  }
  
  auto *mat = dynamic_cast<TPZMatrix<TVar> *>(this->Matrix().operator->());
  if (!mat) {
    std::cerr << "TPZSpoolesSolver::Decompose() - Invalid matrix\n";
    DebugStop();
  }
  
  Decompose(mat);
}

template <class TVar>
void TPZSpoolesSolver<TVar>::Decompose(TPZMatrix<TVar> *mat) {
  auto *spoolesmat = dynamic_cast<TPZFYsmpMatrixSpooles<TVar> *>(mat);
  if (!spoolesmat) {
    std::cerr << "TPZSpoolesSolver::Decompose() - Matrix is not TPZFYsmpMatrixSpooles\n";
    DebugStop();
  }

  FreeSPOOLESMemory();
  fSPOOLESInitialized = true;

  const int neq = static_cast<int>(mat->Rows());
  const int symmetryflag = (fSymmetry == SymProp::Sym) ? 0 : 2; // 0=symmetric, 2=nonsymmetric
  
  if (fMessageLevel > 0) {
    std::cout << "\n========== SPOOLES Factorization ==========\n";
    std::cout << "Matrix size: " << neq << " x " << neq << "\n";
    std::cout << "Symmetry: " << (symmetryflag == 0 ? "Symmetric" : "Non-symmetric") << "\n";
    std::cout << "Threads: " << fNumThreads << "\n";
  }

  // Get COO format from matrix
  TPZVec<int> irn, jcn;
  TPZVec<TVar> val;
  spoolesmat->GetCOOFormat(irn, jcn, val);
  
  const int nnz = irn.size();

  // Create input matrix in COO format
  fInputMatrix = InpMtx_new();
  InpMtx_init(fInputMatrix, INPMTX_BY_ROWS, SPOOLES_REAL, 0, 0);
  
  for (int k = 0; k < nnz; k++) {
    InpMtx_inputRealEntry(fInputMatrix, irn[k], jcn[k], static_cast<double>(val[k]));
  }
  InpMtx_changeStorageMode(fInputMatrix, INPMTX_BY_VECTORS);

  if (fMessageLevel > 2) {
    InpMtx_writeStats(fInputMatrix, stdout);
  }

  // Create graph from the input matrix adjacency
  IVL *adjIVL = InpMtx_fullAdjacency(fInputMatrix);
  int nedges = IVL_tsize(adjIVL);
  
  Graph *graph = Graph_new();
  Graph_init2(graph, 0, neq, 0, nedges, neq, nedges, adjIVL, nullptr, nullptr);
  
  if (fMessageLevel > 2) {
    std::cout << "Graph created with " << neq << " vertices and " << nedges << " edges\n";
  }
  
  // Order the graph using Multiple Minimum Degree
  fETree = orderViaMMD(graph, 12345, fMessageLevel, stdout);
  Graph_free(graph);
  
  if (!fETree || fETree->nfront == 0) {
    std::cerr << "Error: orderViaMMD failed to create elimination tree\n";
    if (fETree) ETree_free(fETree);
    fETree = nullptr;
    DebugStop();
  }
  
  if (fMessageLevel > 1) {
    std::cout << "ETree created with " << fETree->nfront << " fronts, "
              << fETree->nvtx << " vertices\n";
    if (fMessageLevel > 2) {
      ETree_writeStats(fETree, stdout);
    }
  }

  // Create symbolic factorization
  // Get the permutation from the elimination tree
  fOldToNewIV = ETree_oldToNewVtxPerm(fETree);
  int *oldToNew = IV_entries(fOldToNewIV);
  
  fNewToOldIV = ETree_newToOldVtxPerm(fETree);
  
  // Permute the front tree vertices
  ETree_permuteVertices(fETree, fOldToNewIV);
  
  // Permute the input matrix
  InpMtx_permute(fInputMatrix, oldToNew, oldToNew);
  
  // For symmetric matrices, map to upper triangle
  if (symmetryflag != 2) {
    InpMtx_mapToUpperTriangle(fInputMatrix);
  }
  
  // Convert to CHEVRON format required by SymbFac_initFromInpMtx
  InpMtx_changeCoordType(fInputMatrix, INPMTX_BY_CHEVRONS);
  InpMtx_changeStorageMode(fInputMatrix, INPMTX_BY_VECTORS);
  
  fSymbFac = SymbFac_initFromInpMtx(fETree, fInputMatrix);
  
  if (!fSymbFac) {
    std::cerr << "Error: SymbFac_initFromInpMtx failed\n";
    if (fMessageLevel > 0) {
      std::cout << "\nDumping ETree for diagnostics:\n";
      ETree_writeStats(fETree, stdout);
      std::cout << "\nDumping InpMtx for diagnostics:\n";
      InpMtx_writeStats(fInputMatrix, stdout);
    }
    DebugStop();
  }

  // Create front matrix
  fFrontMtx = FrontMtx_new();
  IV *ownersIV = IV_new();
  IV_init(ownersIV, neq, nullptr);
  IV_fill(ownersIV, 0);
  
  fMtxManager = SubMtxManager_new();
  SubMtxManager_init(fMtxManager, NO_LOCK, 0);
  
  FrontMtx_init(fFrontMtx, fETree, fSymbFac, SPOOLES_REAL, symmetryflag, 
                FRONTMTX_DENSE_FRONTS, SPOOLES_PIVOTING, 
                NO_LOCK, 0, ownersIV, fMtxManager,
                fMessageLevel, stdout);
  
  IV_free(ownersIV);

  // Perform factorization
  ChvManager *chvmanager = ChvManager_new();
  ChvManager_init(chvmanager, NO_LOCK, 1);
  
  double tau = 100.0;  // From CalculiX
  double droptol = 0.0;
  int perror = 0;
  double cpus[10];
  int stats[20];
  
  Chv *rootchv = nullptr;
  
  if (fNumThreads == 1) {
    // Serial factorization
    rootchv = FrontMtx_factorInpMtx(fFrontMtx, fInputMatrix, tau, droptol,
                          chvmanager, &perror, cpus, stats,
                          fMessageLevel, stdout);
  } else {
    // Parallel factorization (requires pthreads)
    #ifdef SPOOLES_THREAD_SUPPORT
    rootchv = FrontMtx_MT_factorInpMtx(fFrontMtx, fInputMatrix, tau, droptol,
                             chvmanager, &perror, cpus, stats, fNumThreads, 
                             nullptr, fMessageLevel, stdout);
    #else
    std::cerr << "Warning: SPOOLES compiled without thread support, using serial factorization\n";
    rootchv = FrontMtx_factorInpMtx(fFrontMtx, fInputMatrix, tau, droptol,
                          chvmanager, &perror, cpus, stats,
                          fMessageLevel, stdout);
    #endif
  }

  ChvManager_free(chvmanager);

  if (rootchv != nullptr) {
    std::cerr << "Error: SPOOLES found matrix to be singular\n";
    DebugStop();
  }
  
  if (perror >= 0) {
    std::cerr << "Error in SPOOLES factorization at front " << perror << "\n";
    DebugStop();
  }

  // Post-process the factorization (critical step!)
  FrontMtx_postProcess(fFrontMtx, fMessageLevel, stdout);

  if (fMessageLevel > 0) {
    std::cout << "SPOOLES factorization completed successfully\n";
    std::cout << "===========================================\n\n";
  }

  fDecomposed = true;
}

template <class TVar>
void TPZSpoolesSolver<TVar>::Solve(const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol, 
                                    TPZFMatrix<TVar> *residual) {
  if (!fDecomposed) {
    std::cerr << "TPZSpoolesSolver::Solve() - Matrix not decomposed\n";
    Decompose();
  }

  auto *mat = dynamic_cast<TPZMatrix<TVar> *>(this->Matrix().operator->());
  if (!mat) {
    std::cerr << "TPZSpoolesSolver::Solve() - Invalid matrix\n";
    DebugStop();
  }

  Solve(mat, rhs, sol);

  if (residual) {
    this->Matrix()->Residual(sol, rhs, *residual);
  }
}

template <class TVar>
void TPZSpoolesSolver<TVar>::Solve(const TPZMatrix<TVar> *mat, const TPZFMatrix<TVar> &rhs, 
                                    TPZFMatrix<TVar> &sol) const {
  const int neq = static_cast<int>(mat->Rows());
  const int nrhs = rhs.Cols();

  sol.Resize(neq, nrhs);

  // Create RHS matrix
  DenseMtx *mtxY = DenseMtx_new();
  DenseMtx_init(mtxY, SPOOLES_REAL, 0, 0, neq, nrhs, 1, neq);
  
  // Copy RHS data
  double *entriesY = DenseMtx_entries(mtxY);
  for (int j = 0; j < nrhs; j++) {
    for (int i = 0; i < neq; i++) {
      entriesY[i + j*neq] = static_cast<double>(rhs.GetVal(i, j));
    }
  }

  // std::cout << "DEBUG: RHS before manual permutation:\n";
  // for (int i = 0; i < std::min(neq, 9); i++) {
  //   std::cout << "  rhs[" << i << "] = " << entriesY[i] << "\n";
  // }

  // Manual permutation: newRHS[oldToNew[i]] = oldRHS[i]
  TPZVec<double> tempRHS(neq);
  int *oldToNew = IV_entries(fOldToNewIV);
  // std::cout << "DEBUG: Permutation mapping:\n";
  // for (int i = 0; i < neq; i++) {
  //   tempRHS[oldToNew[i]] = entriesY[i];
  //   std::cout << "  tempRHS[" << oldToNew[i] << "] = oldRHS[" << i << "] = " << entriesY[i] << "\n";
  // }
  for (int i = 0; i < neq; i++) {
    entriesY[i] = tempRHS[i];
  }
  
  // std::cout << "DEBUG: RHS after manual permutation:\n";
  // for (int i = 0; i < std::min(neq, 9); i++) {
  //   std::cout << "  rhs[" << i << "] = " << entriesY[i] << "\n";
  // }

  // Create solution matrix
  DenseMtx *mtxX = DenseMtx_new();
  DenseMtx_init(mtxX, SPOOLES_REAL, 0, 0, neq, nrhs, 1, neq);
  DenseMtx_zero(mtxX);

  // Solve
  double cpus[10];
  FrontMtx_solve(fFrontMtx, mtxX, mtxY, fMtxManager, cpus,
                 fMessageLevel, stdout);

  double *entriesXbeforePerm = DenseMtx_entries(mtxX);
  // std::cout << "DEBUG: Solution before manual permutation back:\n";
  // for (int i = 0; i < std::min(neq, 9); i++) {
  //   std::cout << "  sol[" << i << "] = " << entriesXbeforePerm[i] << "\n";
  // }

  // Manual permutation back: newSol[i] = oldSol[newToOld[i]]
  TPZVec<double> tempSol(neq);
  int *newToOld = IV_entries(fNewToOldIV);
  for (int i = 0; i < neq; i++) {
    tempSol[i] = entriesXbeforePerm[newToOld[i]];
  }
  for (int i = 0; i < neq; i++) {
    entriesXbeforePerm[i] = tempSol[i];
  }
  
  // std::cout << "DEBUG: Solution after manual permutation back:\n";
  // for (int i = 0; i < std::min(neq, 9); i++) {
  //   std::cout << "  sol[" << i << "] = " << entriesXbeforePerm[i] << "\n";
  // }

  // Copy solution back (get entries AFTER permutation!)
  double *entriesX = DenseMtx_entries(mtxX);
  for (int j = 0; j < nrhs; j++) {
    for (int i = 0; i < neq; i++) {
      sol.PutVal(i, j, static_cast<TVar>(entriesX[i + j*neq]));
    }
  }

  DenseMtx_free(mtxY);
  DenseMtx_free(mtxX);
}

// Explicit instantiations (SPOOLES only supports real types)
template class TPZSpoolesSolver<float>;
template class TPZSpoolesSolver<double>;
template class TPZSpoolesSolver<long double>;
