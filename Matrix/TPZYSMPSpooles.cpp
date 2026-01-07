/**
 * @file
 * @brief Contains the implementation of the TPZFYsmpMatrixSpooles methods.
 */

#include "TPZYSMPSpooles.h"
#include "pzfmatrix.h"

// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZFYsmpMatrixSpooles<TVar>::TPZFYsmpMatrixSpooles() : 
    TPZRegisterClassId(&TPZFYsmpMatrixSpooles::ClassId),
    TPZFYsmpMatrix<TVar>() {
}

template<class TVar>
TPZFYsmpMatrixSpooles<TVar>::TPZFYsmpMatrixSpooles(const TPZFYsmpMatrixSpooles<TVar> &cp) :
    TPZRegisterClassId(&TPZFYsmpMatrixSpooles::ClassId),
    TPZFYsmpMatrix<TVar>(cp) {
    // Note: fSpoolesControl is NOT copied - each matrix gets a fresh SPOOLES instance
}

template<class TVar>
TPZFYsmpMatrixSpooles<TVar>::TPZFYsmpMatrixSpooles(TPZFYsmpMatrixSpooles<TVar> &&cp) :
    TPZRegisterClassId(&TPZFYsmpMatrixSpooles::ClassId),
    TPZFYsmpMatrix<TVar>(std::move(cp)),
    fSpoolesControl(std::move(cp.fSpoolesControl)) {
}

template<class TVar>
TPZFYsmpMatrixSpooles<TVar>& TPZFYsmpMatrixSpooles<TVar>::operator=(const TPZFYsmpMatrixSpooles<TVar> &copy) {
    if (this != &copy) {
        TPZFYsmpMatrix<TVar>::operator=(copy);
        // Note: fSpoolesControl is NOT copied - each matrix gets a fresh SPOOLES instance
        fSpoolesControl = TPZSpoolesSolver<TVar>();
        // Reset decomposition state
        this->SetIsDecomposed(ENoDecompose);
    }
    return *this;
}

template<class TVar>
TPZFYsmpMatrixSpooles<TVar>& TPZFYsmpMatrixSpooles<TVar>::operator=(TPZFYsmpMatrixSpooles<TVar> &&copy) {
    if (this != &copy) {
        TPZFYsmpMatrix<TVar>::operator=(std::move(copy));
        fSpoolesControl = std::move(copy.fSpoolesControl);
    }
    return *this;
}

template<class TVar>
void TPZFYsmpMatrixSpooles<TVar>::CopyFrom(const TPZMatrix<TVar> *mat) {
    auto *from = dynamic_cast<const TPZFYsmpMatrixSpooles<TVar> *>(mat);
    if (from) {
        *this = *from;
        return;
    }
    
    auto *fromBase = dynamic_cast<const TPZFYsmpMatrix<TVar> *>(mat);
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
int TPZFYsmpMatrixSpooles<TVar>::ClassId() const {
    return Hash("TPZFYsmpMatrixSpooles") ^ TPZFYsmpMatrix<TVar>::ClassId() << 1;
}

template<class TVar>
void TPZFYsmpMatrixSpooles<TVar>::MultAdd(const TPZFMatrix<TVar> &x,
                                          const TPZFMatrix<TVar> &y,
                                          TPZFMatrix<TVar> &z,
                                          const TVar alpha, const TVar beta,
                                          const int opt) const {
    // Use parent class implementation (standard CSR matrix-vector multiplication)
    TPZFYsmpMatrix<TVar>::MultAdd(x, y, z, alpha, beta, opt);
}

template<class TVar>
void TPZFYsmpMatrixSpooles<TVar>::SetIsDecomposed(DecomposeType val) {
    this->fDecomposed = val;
    if (val == ENoDecompose) {
        fSpoolesControl.FreeSPOOLESMemory();
    }
}

template<class TVar>
int TPZFYsmpMatrixSpooles<TVar>::Decompose(const DecomposeType dt) {
    if (this->fDecomposed != ENoDecompose) {
        return 1;
    }

    // Set matrix type based on decomposition type
    if (dt == ELDLt) {
        fSpoolesControl.SetMatrixType(SymProp::Sym, 
                                       TPZSpoolesSolver<TVar>::MProperty::EIndefinite);
    } else if (dt == ECholesky) {
        fSpoolesControl.SetMatrixType(SymProp::Sym, 
                                       TPZSpoolesSolver<TVar>::MProperty::EPositiveDefinite);
    } else if (dt == ELU) {
        fSpoolesControl.SetMatrixType(SymProp::NonSym, 
                                       TPZSpoolesSolver<TVar>::MProperty::EIndefinite);
    } else {
        std::cerr << "TPZFYsmpMatrixSpooles::Decompose - Unsupported decomposition type\n";
        DebugStop();
    }

    fSpoolesControl.SetMatrix(this);
    fSpoolesControl.Decompose();
    
    this->fDecomposed = dt;
    return 1;
}

template<class TVar>
int TPZFYsmpMatrixSpooles<TVar>::SolveDirect(TPZFMatrix<TVar>& F, const DecomposeType dt) {
    if (this->fDecomposed != dt) {
        Decompose(dt);
    }
    
    fSpoolesControl.Solve(F, F);
    return 1;
}

template<class TVar>
int TPZFYsmpMatrixSpooles<TVar>::SolveDirect(TPZFMatrix<TVar>& F, const DecomposeType dt) const {
    if (this->fDecomposed != dt) {
        std::cerr << "TPZFYsmpMatrixSpooles::SolveDirect (const) - Matrix not decomposed\n";
        DebugStop();
    }
    
    // Cast away const for solve (SPOOLES modifies internal state)
    auto *non_const_this = const_cast<TPZFYsmpMatrixSpooles<TVar>*>(this);
    non_const_this->fSpoolesControl.Solve(F, F);
    return 1;
}

template<class TVar>
void TPZFYsmpMatrixSpooles<TVar>::GetCOOFormat(TPZVec<int> &irn, TPZVec<int> &jcn, 
                                                TPZVec<TVar> &val) const {
    // Convert CSR to COO format (0-based indexing for SPOOLES)
    const int64_t nrows = this->Rows();
    const auto &ia = this->fIA;
    const auto &ja = this->fJA;
    const auto &a = this->fA;
    
    // Count non-zeros
    int64_t nnz = ia[nrows];
    
    irn.Resize(nnz);
    jcn.Resize(nnz);
    val.Resize(nnz);
    
    int64_t k = 0;
    for (int64_t i = 0; i < nrows; i++) {
        for (int64_t j = ia[i]; j < ia[i+1]; j++) {
            irn[k] = static_cast<int>(i);
            jcn[k] = static_cast<int>(ja[j]);
            val[k] = a[j];
            k++;
        }
    }
}

// Explicit instantiations (SPOOLES only supports real types)
template class TPZFYsmpMatrixSpooles<float>;
template class TPZFYsmpMatrixSpooles<double>;
template class TPZFYsmpMatrixSpooles<long double>;
