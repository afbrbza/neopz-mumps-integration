/**
 * @file
 * @brief Contains TPZFYsmpMatrixSpooles class for SPOOLES sparse matrix format
 */

#ifndef YSMPMAT_SPOOLES_H
#define YSMPMAT_SPOOLES_H

#include "TPZYSMPMatrix.h"
#include "TPZSpoolesSolver.h"

/**
 * @brief Implements a sparse matrix using SPOOLES solver. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZFYsmpMatrixSpooles : public TPZFYsmpMatrix<TVar>{
    
    friend class TPZSpoolesSolver<TVar>;
    
public:

  TPZFYsmpMatrixSpooles();
  
  /** @brief Constructors from parent class*/
  using TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix;
    
  /** @brief Copy constructor - creates new SPOOLES instance */
  TPZFYsmpMatrixSpooles(const TPZFYsmpMatrixSpooles<TVar> &cp);
  
  /** @brief Move constructor*/
  TPZFYsmpMatrixSpooles(TPZFYsmpMatrixSpooles<TVar> &&cp);
  
  /** @brief Copy-assignment operator - creates new SPOOLES instance*/
  TPZFYsmpMatrixSpooles &operator=(const TPZFYsmpMatrixSpooles<TVar> &copy);
  
  /** @brief Move-assignment operator*/
  TPZFYsmpMatrixSpooles &operator=(TPZFYsmpMatrixSpooles<TVar> &&copy);

  /** @brief Copy constructor from generic sparse matrix*/
  TPZFYsmpMatrixSpooles(const TPZFYsmpMatrix<TVar> &cp)
    : TPZFYsmpMatrix<TVar>(cp) {}
    
  /** @brief Move constructor from generic sparse matrix*/
  TPZFYsmpMatrixSpooles(TPZFYsmpMatrix<TVar> &&rval)
    : TPZFYsmpMatrix<TVar>(rval) {}
    
  /** @brief Copy-assignment operator from generic sparse matrix*/
  TPZFYsmpMatrixSpooles &operator=(const TPZFYsmpMatrix<TVar> &cp)
  { TPZFYsmpMatrix<TVar>::operator=(cp); return *this;}
  
  /** @brief Move-assignment operator from generic sparse matrix*/
  TPZFYsmpMatrixSpooles &operator=(TPZFYsmpMatrix<TVar> &&rval)
  { TPZFYsmpMatrix<TVar>::operator=(rval); return *this;}
  
  inline TPZFYsmpMatrixSpooles<TVar>*NewMatrix() const override {
    return new TPZFYsmpMatrixSpooles<TVar>{};
  }
  
  CLONEDEF(TPZFYsmpMatrixSpooles)
  
  /** @brief Destructor */
  ~TPZFYsmpMatrixSpooles() = default;

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
  
  //! Get COO format arrays (0-based indexing)
  void GetCOOFormat(TPZVec<int> &irn, TPZVec<int> &jcn, TPZVec<TVar> &val) const;
  
private:
  
  TPZSpoolesSolver<TVar> fSpoolesControl;
};

#endif // YSMPMAT_SPOOLES_H
