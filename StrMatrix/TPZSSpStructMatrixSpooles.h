/**
 * @file
 * @brief Contains the TPZSSpStructMatrixSpooles class which implements sparse symmetric structural matrices.
 */

#ifndef TPZSSpStructMatrixSpooles_H
#define TPZSSpStructMatrixSpooles_H

#include "pzstack.h"
#include "pzstrmatrixor.h"
#include "TPZSSpStructMatrix.h"

/**
 * @brief Implements a sparse symmetric structural matrix using TPZSYsmpMatrixSpooles as a storage format.
 * This specialized version uses SPOOLES solver.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSSpStructMatrixSpooles : public TPZSSpStructMatrix<TVar,TPar> {
public:
    // Inherit constructors from parent class
    using TPZSSpStructMatrix<TVar,TPar>::TPZSSpStructMatrix;
    
protected:
    // Override SetupMatrixData if you need custom SPOOLES behavior
    virtual TPZMatrix<TVar> * Create() override;
    virtual TPZSSpStructMatrixSpooles * Clone() override;
    virtual TPZMatrix<TVar> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex) override;
    friend TPZPersistenceManager;
};

#endif //TPZSSpStructMatrixSpooles_H
