/**
 * @file
 * @brief Contains the TPZSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSPSTRUCTMATRIXMUMPS_H
#define TPZSPSTRUCTMATRIXMUMPS_H

#include "TPZSpStructMatrix.h"
#include "pzstack.h"
#include "pzstrmatrixor.h"
/**
 * @brief Implements a sparse structural matrix using TPZFYsmpMatrix as a storage format.
 * @ingroup structural
 */
template <class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSpStructMatrixMumps : public TPZSpStructMatrix<TVar, TPar> {
public:
  using TPZSpStructMatrix<TVar, TPar>::TPZSpStructMatrix;

protected:
  virtual TPZMatrix<TVar> *SetupMatrixData(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex);
  friend TPZPersistenceManager;
};

#endif // TPZSPSTRUCTMATRIXMUMPS_H
