/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrixSpooles methods.
 */

#include "TPZSYSMPSpooles.h"
#include "pzfmatrix.h"

// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZSYsmpMatrixSpooles<TVar>::TPZSYsmpMatrixSpooles() : 
    TPZRegisterClassId(&TPZSYsmpMatrixSpooles::ClassId),
    TPZSYsmpMatrix<TVar>() {
}

template<class TVar>
TPZSYsmpMatrixSpooles<TVar>::TPZSYsmpMatrixSpooles(const TPZSYsmpMatrixSpooles<TVar> &cp) :
    TPZRegisterClassId(&TPZSYsmpMatrixSpooles::ClassId),
    TPZSYsmpMatrix<TVar>(cp) {
    // Note: fSpoolesControl is NOT copied - each matrix gets a fresh SPOOLES instance
}

template<class TVar>
TPZSYsmpMatrixSpooles<TVar>::TPZSYsmpMatrixSpooles(TPZSYsmpMatrixSpooles<TVar> &&cp) :
    TPZRegisterClassId(&TPZSYsmpMatrixSpooles::ClassId),
    TPZSYsmpMatrix<TVar>(std::move(cp)),
    fSpoolesControl(std::move(cp.fSpoolesControl)) {
}

template<class TVar>
TPZSYsmpMatrixSpooles<TVar>& TPZSYsmpMatrixSpooles<TVar>::operator=(const TPZSYsmpMatrixSpooles<TVar> &copy) {
    if (this != &copy) {
        TPZSYsmpMatrix<TVar>::operator=(copy);
        // Note: fSpoolesControl is NOT copied - each matrix gets a fresh SPOOLES instance
        fSpoolesControl = TPZSpoolesSolver<TVar>();
        // Reset decomposition state
        this->SetIsDecomposed(ENoDecompose);
    }
    return *this;
}

template<class TVar>
TPZSYsmpMatrixSpooles<TVar>& TPZSYsmpMatrixSpooles<TVar>::operator=(TPZSYsmpMatrixSpooles<TVar> &&copy) {
    if (this != &copy) {
        TPZSYsmpMatrix<TVar>::operator=(std::move(copy));
        fSpoolesControl = std::move(copy.fSpoolesControl);
    }
    return *this;
}

template<class TVar>
void TPZSYsmpMatrixSpooles<TVar>::CopyFrom(const TPZMatrix<TVar> *mat) {
    auto *from = dynamic_cast<const TPZSYsmpMatrixSpooles<TVar> *>(mat);
    if (from) {
        *this = *from;
        return;
    }
    
    auto *fromBase = dynamic_cast<const TPZSYsmpMatrix<TVar> *>(mat);
    if (fromBase && fromBase->IsDecomposed() == ENoDecompose) {
        *this = *fromBase;
        return;
    }
    
    PZError << __PRETTY_FUNCTION__;
    PZError << "\nERROR: Called with incompatible type\n";
    PZError << "Aborting...\n";
    DebugStop();
}

template<class TVar>
int TPZSYsmpMatrixSpooles<TVar>::ClassId() const {
    return Hash("TPZSYsmpMatrixSpooles") ^ TPZSYsmpMatrix<TVar>::ClassId() << 1;
}

template<class TVar>
void TPZSYsmpMatrixSpooles<TVar>::MultAdd(const TPZFMatrix<TVar> &x,
                                          const TPZFMatrix<TVar> &y,
                                          TPZFMatrix<TVar> &z,
                                          const TVar alpha, const TVar beta,
                                          const int opt) const {
    // Use parent class implementation (standard symmetric sparse matrix-vector multiplication)
    TPZSYsmpMatrix<TVar>::MultAdd(x, y, z, alpha, beta, opt);
}

template<class TVar>
void TPZSYsmpMatrixSpooles<TVar>::SetIsDecomposed(DecomposeType val) {
    this->fDecomposed = val;
    if (val == ENoDecompose) {
        fSpoolesControl.FreeSPOOLESMemory();
    }
}

template<class TVar>
int TPZSYsmpMatrixSpooles<TVar>::Decompose(const DecomposeType dt) {
    if (this->fDecomposed != ENoDecompose) {
        return 1;
    }

    // Set matrix type based on decomposition type
    // Symmetric matrices use SymProp::Sym
    if (dt == ELDLt) {
        fSpoolesControl.SetMatrixType(SymProp::Sym, 
                                       TPZSpoolesSolver<TVar>::MProperty::EIndefinite);
    } else if (dt == ECholesky) {
        fSpoolesControl.SetMatrixType(SymProp::Sym, 
                                       TPZSpoolesSolver<TVar>::MProperty::EPositiveDefinite);
    } else {
        std::cerr << "TPZSYsmpMatrixSpooles::Decompose - Unsupported decomposition type for symmetric matrix\n";
        std::cerr << "Use ECholesky or ELDLt\n";
        DebugStop();
    }

    fSpoolesControl.SetMatrix(this);
    fSpoolesControl.Decompose();
    
    this->fDecomposed = dt;
    return 1;
}

template<class TVar>
int TPZSYsmpMatrixSpooles<TVar>::SolveDirect(TPZFMatrix<TVar>& F, const DecomposeType dt) {
    if (this->fDecomposed != dt) {
        Decompose(dt);
    }
    
    fSpoolesControl.Solve(F, F);
    return 1;
}

template<class TVar>
int TPZSYsmpMatrixSpooles<TVar>::SolveDirect(TPZFMatrix<TVar>& F, const DecomposeType dt) const {
    if (this->fDecomposed != dt) {
        std::cerr << "TPZSYsmpMatrixSpooles::SolveDirect (const) - Matrix not decomposed\n";
        DebugStop();
    }
    
    // Cast away const for solve (SPOOLES modifies internal state)
    auto *non_const_this = const_cast<TPZSYsmpMatrixSpooles<TVar>*>(this);
    non_const_this->fSpoolesControl.Solve(F, F);
    return 1;
}

template<class TVar>
void TPZSYsmpMatrixSpooles<TVar>::GetCOOFormat(TPZVec<int> &irn, TPZVec<int> &jcn, 
                                                TPZVec<TVar> &val) const {
    // Convert symmetric CSR to COO format (upper triangular only, 0-based indexing for SPOOLES)
    const int64_t nrows = this->Rows();
    const auto &ia = this->fIA;
    const auto &ja = this->fJA;
    const auto &a = this->fA;
    
    // Count non-zeros in upper triangular part
    int64_t nnz = 0;
    for (int64_t i = 0; i < nrows; i++) {
        for (int64_t j = ia[i]; j < ia[i+1]; j++) {
            if (ja[j] >= i) {  // Only upper triangular (including diagonal)
                nnz++;
            }
        }
    }
    
    irn.Resize(nnz);
    jcn.Resize(nnz);
    val.Resize(nnz);
    
    int64_t k = 0;
    for (int64_t i = 0; i < nrows; i++) {
        for (int64_t j = ia[i]; j < ia[i+1]; j++) {
            if (ja[j] >= i) {  // Only upper triangular
                irn[k] = static_cast<int>(i);
                jcn[k] = static_cast<int>(ja[j]);
                val[k] = a[j];
                k++;
            }
        }
    }
}

// Explicit instantiations (SPOOLES only supports real types)
template class TPZSYsmpMatrixSpooles<float>;
template class TPZSYsmpMatrixSpooles<double>;
template class TPZSYsmpMatrixSpooles<long double>;
