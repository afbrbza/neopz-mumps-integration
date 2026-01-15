/**
 * @file
 * @brief Contains the TPZSSpStructMatrixMumps class which implements sparse structural matrices.
 */

#ifndef TPZSSpStructMatrixMumps_H
#define TPZSSpStructMatrixMumps_H

#include "pzstack.h"

#include "pzstrmatrixor.h"
#include "TPZSSpStructMatrix.h"
/**
 * @brief Implements a sparse symmetric structural matrix using TPZSYsmpMatrix as a storage format.
 * This specialized version uses MUMPS solver.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSSpStructMatrixMumps : public TPZSSpStructMatrix<TVar,TPar> {
public:
    // Inherit constructors from parent class
    using TPZSSpStructMatrix<TVar,TPar>::TPZSSpStructMatrix;
    
protected:
    // Override SetupMatrixData if you need custom MUMPS behavior
    virtual TPZMatrix<TVar> * Create() override;
	virtual TPZSSpStructMatrixMumps * Clone() override;
    virtual TPZMatrix<TVar> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex) override;
};

#endif //TPZSSpStructMatrixMumps_H
