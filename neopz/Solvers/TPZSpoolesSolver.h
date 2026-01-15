#ifndef TPZSPOOLEESSOLVER_H
#define TPZSPOOLEESSOLVER_H

#include "TPZMatrixSolver.h"
#include "pzmanvector.h"
#include "tpzautopointer.h"

// SPOOLES headers
extern "C" {
#include "misc.h"
#include "InpMtx.h"
#include "ETree.h"
#include "IVL.h"
#include "IV.h"
#include "Graph.h"
#include "FrontMtx.h"
#include "ChvManager.h"
#include "SubMtxManager.h"
#include "DenseMtx.h"
#include "SymbFac.h"
}

template <class TVar>
class TPZFYsmpMatrixSpooles;

template <typename TVar>
class TPZSpoolesSolver : public TPZMatrixSolver<TVar> {

  friend class TPZFYsmpMatrixSpooles<TVar>;

public:
  enum class MProperty { 
    ENonInitialized = 0,
    EPositiveDefinite,
    EIndefinite 
  };

  TPZSpoolesSolver();

  // Delete copy operations - SPOOLES data structures cannot be safely copied
  TPZSpoolesSolver(const TPZSpoolesSolver &copy) = delete;
  TPZSpoolesSolver &operator=(const TPZSpoolesSolver &copy) = delete;

  // Move operations are allowed
  TPZSpoolesSolver(TPZSpoolesSolver &&copy) noexcept;
  TPZSpoolesSolver &operator=(TPZSpoolesSolver &&copy) noexcept;

  virtual ~TPZSpoolesSolver();

  void SetMatrix(TPZAutoPointer<TPZBaseMatrix> Refmat) override;

  void Decompose() override;

  void Solve(const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol, TPZFMatrix<TVar> *residual = nullptr) override;

  void SetMatrixType(SymProp symtype, MProperty prop);

  void SetMessageLevel(int lvl);

  int GetMessageLevel() const { return fMessageLevel; }

  TPZSpoolesSolver<TVar> *Clone() const override;

  void SetNumThreads(int nthreads) { fNumThreads = nthreads; }
  
  int GetNumThreads() const { return fNumThreads; }

  void FreeSPOOLESMemory();

protected:
  void Decompose(TPZMatrix<TVar> *mat);

  void Solve(const TPZMatrix<TVar> *mat, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const;

  SymProp fSymmetry{SymProp::NonSym};

  MProperty fProperty{MProperty::ENonInitialized};

  // SPOOLES data structures
  InpMtx *fInputMatrix{nullptr};
  ETree *fETree{nullptr};
  IVL *fSymbFac{nullptr};
  FrontMtx *fFrontMtx{nullptr};
  SubMtxManager *fMtxManager{nullptr};
  DenseMtx *fSolutionMtx{nullptr};
  DenseMtx *fRhsMtx{nullptr};
  IV *fOldToNewIV{nullptr};
  IV *fNewToOldIV{nullptr};

  int fMessageLevel{0};
  int fNumThreads{1};
  
  int fMatrixType{0}; // 0=real, 1=complex

  bool fDecomposed{false};

  bool fSPOOLESInitialized{false};
};

#endif /* TPZSPOOLEESSOLVER_H */
