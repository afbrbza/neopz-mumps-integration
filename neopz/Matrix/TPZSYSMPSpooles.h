/**
 * @file
 * @brief Contains TPZSYsmpMatrixSpooles class for symmetric sparse matrices with SPOOLES
 */

#ifndef SYSMPMAT_SPOOLES_H
#define SYSMPMAT_SPOOLES_H

#include "TPZSYSMPMatrix.h"
#include "TPZSpoolesSolver.h"

/**
 * @brief Implements a symmetric sparse matrix using SPOOLES solver. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSYsmpMatrixSpooles : public TPZSYsmpMatrix<TVar>{
    
    friend class TPZSpoolesSolver<TVar>;
    
public:

  TPZSYsmpMatrixSpooles();
  
  /** @brief Constructors from parent class*/
  using TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix;
    
  /** @brief Copy constructor - creates new SPOOLES instance */
  TPZSYsmpMatrixSpooles(const TPZSYsmpMatrixSpooles<TVar> &cp);
  
  /** @brief Move constructor*/
  TPZSYsmpMatrixSpooles(TPZSYsmpMatrixSpooles<TVar> &&cp);
  
  /** @brief Copy-assignment operator - creates new SPOOLES instance*/
  TPZSYsmpMatrixSpooles &operator=(const TPZSYsmpMatrixSpooles<TVar> &copy);
  
  /** @brief Move-assignment operator*/
  TPZSYsmpMatrixSpooles &operator=(TPZSYsmpMatrixSpooles<TVar> &&copy);

  /** @brief Copy constructor from generic sparse matrix*/
  TPZSYsmpMatrixSpooles(const TPZSYsmpMatrix<TVar> &cp)
    : TPZSYsmpMatrix<TVar>(cp) {}
    
  /** @brief Move constructor from generic sparse matrix*/
  TPZSYsmpMatrixSpooles(TPZSYsmpMatrix<TVar> &&rval)
    : TPZSYsmpMatrix<TVar>(rval) {}
    
  /** @brief Copy-assignment operator from generic sparse matrix*/
  TPZSYsmpMatrixSpooles &operator=(const TPZSYsmpMatrix<TVar> &cp)
  { TPZSYsmpMatrix<TVar>::operator=(cp); return *this;}
  
  /** @brief Move-assignment operator from generic sparse matrix*/
  TPZSYsmpMatrixSpooles &operator=(TPZSYsmpMatrix<TVar> &&rval)
  { TPZSYsmpMatrix<TVar>::operator=(rval); return *this;}
  
  inline TPZSYsmpMatrixSpooles<TVar>*NewMatrix() const override {
    return new TPZSYsmpMatrixSpooles<TVar>{};
  }
  
  CLONEDEF(TPZSYsmpMatrixSpooles)
  
  /** @brief Destructor */
  ~TPZSYsmpMatrixSpooles() = default;

  /** @brief Creates a copy from another sparse matrix*/
  void CopyFrom(const TPZMatrix<TVar> *mat) override;

  void MultAdd(const TPZFMatrix<TVar> &x, const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
               const TVar alpha=1., const TVar beta = 0., const int opt = 0) const override;
  
  void SetIsDecomposed(DecomposeType val) override;
  
  virtual int Decompose(const DecomposeType dt) override;
  
  virtual int SolveDirect(TPZFMatrix<TVar>& F, const DecomposeType dt) override;
  
  virtual int SolveDirect(TPZFMatrix<TVar>& F, const DecomposeType dt) const override;

  int ClassId() const override;

  //! Gets reference to TPZSpoolesSolver instance for fine-tuning
  TPZSpoolesSolver<TVar> & GetSpoolesControl() {
    return fSpoolesControl;
  }
  
  //! Get COO format arrays for symmetric matrix (upper triangular only, 0-based indexing)
  void GetCOOFormat(TPZVec<int> &irn, TPZVec<int> &jcn, TPZVec<TVar> &val) const;
  
private:
  
  TPZSpoolesSolver<TVar> fSpoolesControl;
};

#endif // SYSMPMAT_SPOOLES_H
