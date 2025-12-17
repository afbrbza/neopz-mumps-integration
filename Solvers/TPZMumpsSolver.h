#ifndef TPZMUMPSSOLVER_H
#define TPZMUMPSSOLVER_H

#include "TPZMatrixSolver.h"

#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "dmumps_c.h"

#define JOB_INIT -1
#define JOB_END -2
#define JOB_ANALYSIS 1
#define JOB_FACTORIZE 2
#define JOB_SOLVE 3

#define USE_COMM_WORLD -987654 // sequential run without MPI
#define MUMPS_HOST_PAR 1

template <class TVar>
class TPZFYsmpMatrixMumps;

template <class TVar>
class TPZSYsmpMatrixMumps;

template <typename TVar>
class TPZMumpsSolver : public TPZMatrixSolver<TVar> {

  // they need access to SetRawMatrix and ReallySolve
  friend class TPZFYsmpMatrixMumps<TVar>;
  friend class TPZSYsmpMatrixMumps<TVar>;

public:
  enum class MProperty { ENonInitialized = 0,
                         EPositiveDefinite,
                         EIndefinite };

  TPZMumpsSolver();

  // Delete copy operations - MUMPS data structures cannot be safely copied
  TPZMumpsSolver(const TPZMumpsSolver &copy) = delete;
  TPZMumpsSolver &operator=(const TPZMumpsSolver &copy) = delete;

  // Move operations are allowed
  TPZMumpsSolver(TPZMumpsSolver &&copy) noexcept;
  TPZMumpsSolver &operator=(TPZMumpsSolver &&copy) noexcept;

  virtual ~TPZMumpsSolver();

  void SetMatrix(TPZAutoPointer<TPZBaseMatrix> Refmat) override;

  void Decompose() override;

  void Solve(const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol, TPZFMatrix<TVar> *residual = nullptr) override;

  void SetMatrixType(SymProp symtype, MProperty prop);

  void SetMatrixType(long long type) { fMatrixType = type; }

  void SetMessageLevel(int lvl);

  int GetMessageLevel() const { return fMessageLevel; }

  TPZMumpsSolver<TVar> *Clone() const override;

  [[nodiscard]] inline TPZVec<long long> GetParam() const { return fParam; }

  void SetParam(const TPZVec<long long> &p);

  void ResetParam();

  [[nodiscard]] bool HasCustomSettings() const { return fCustomSettings; }

  TPZVec<long long> &GetPermutationVec() { return fPermutation; }

  void FreeMumpsMemory();

protected:
  long long MatrixType();

  void Decompose(TPZMatrix<TVar> *mat);

  void Solve(const TPZMatrix<TVar> *mat, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const;

  SymProp fSymmetry{SymProp::NonSym};

  MProperty fProperty{MProperty::ENonInitialized};

  mutable DMUMPS_STRUC_C fMumpsData;

  TPZManVector<long long, 64> fParam;

  // Persistent storage for 1-based index arrays expected by MUMPS
  TPZManVector<int64_t> fIRN1Based;
  TPZManVector<int64_t> fJCN1Based;

  long long fMax_num_factors{1};

  long long fMatrix_num{1};

  long long fMessageLevel{0};

  long long fError{0};

  TPZVec<long long> fPermutation;

  long long fMatrixType{0};

  bool fDecomposed{false};

  bool fMumpsInitialized{false};

  bool fCustomSettings{false};
};

#endif /* TPZMUMPSSOLVER_H */
