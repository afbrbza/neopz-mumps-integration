/**
 * @file
 * @brief Contains the TPZSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSPSTRUCTMATRIX_SPOOLES_H
#define TPZSPSTRUCTMATRIX_SPOOLES_H

#include "TPZSpStructMatrix.h"
#include "pzstack.h"
#include "pzstrmatrixor.h"

/**
 * @brief Implements a sparse structural matrix using TPZFYsmpMatrixSpooles as a storage format.
 * @ingroup structural
 */
template <class TVar = STATE, class TPar = TPZStructMatrixOR<TVar>>
class TPZSpStructMatrixSpooles : public TPZSpStructMatrix<TVar, TPar> {
public:
  using TPZSpStructMatrix<TVar, TPar>::TPZSpStructMatrix;

protected:
  virtual TPZMatrix<TVar> *Create() override;
  virtual TPZSpStructMatrixSpooles *Clone() override;
  virtual TPZMatrix<TVar> *SetupMatrixData(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex);
  friend TPZPersistenceManager;
};

#endif // TPZSPSTRUCTMATRIX_SPOOLES_H
