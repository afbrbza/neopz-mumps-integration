/**
 * @file
 * @brief Contains the implementation of the TPZSSpStructMatrixSpooles methods. 
 */

#include "TPZSSpStructMatrixSpooles.h"
#include "pzcmesh.h"
#include "TPZSYSMPMatrix.h"
#include "TPZRenumbering.h"
#include "TPZSYSMPSpooles.h"
#include "TPZTimer.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar, class TPar>
TPZSSpStructMatrixSpooles<TVar, TPar> * TPZSSpStructMatrixSpooles<TVar,TPar>::Clone(){
    return new TPZSSpStructMatrixSpooles(*this);
}

template <class TVar, class TPar>
TPZMatrix<TVar> *TPZSSpStructMatrixSpooles<TVar, TPar>::Create() {
  int64_t neq = this->fMesh->NEquations();

  TPZStack<int64_t> elgraph;
  TPZVec<int64_t> elgraphindex;
  this->fMesh->ComputeElGraph(elgraph, elgraphindex);

  TPZMatrix<TVar> *mat = SetupMatrixData(elgraph, elgraphindex);
  return mat;
}

template<class TVar, class TPar>
TPZMatrix<TVar> *
TPZSSpStructMatrixSpooles<TVar,TPar>::SetupMatrixData(TPZStack<int64_t> & elgraph,
                                              TPZVec<int64_t> &elgraphindex){
    
    const int64_t neq = this->fEquationFilter.NActiveEquations();
    TPZSYsmpMatrixSpooles<TVar> * mat = new TPZSYsmpMatrixSpooles<TVar>(neq,neq);
    
    /**Creates a element graph*/
    TPZRenumbering graph;
    graph.SetElementsNodes(elgraphindex.NElements() -1 ,
                           this->fMesh->NIndependentConnects());
	
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    graph.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    const int64_t nblock = nodegraphindex.NElements()-1;
    // number of values in the sparse matrix (symmetric, upper triangular only)
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(int64_t i=0;i<nblock;i++){
		const int64_t iblsize = this->fMesh->Block().Size(i);
		const int64_t iblpos = this->fMesh->Block().Position(i);
        const int64_t numactive = this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
		totaleq += iblsize;
		// Diagonal block
		totalvar += (iblsize * (iblsize + 1)) / 2;  // Only upper triangular
		
		const int64_t icfirst = nodegraphindex[i];
		const int64_t iclast = nodegraphindex[i+1];
		for(auto j=icfirst;j<iclast;j++) {
			const int64_t col = nodegraph[j];
			if (col > i) {  // Only upper triangular
				const int64_t colsize = this->fMesh->Block().Size(col);
				const int64_t colpos = this->fMesh->Block().Position(col);
				const int64_t numactive = this->fEquationFilter.NumActive(colpos, colpos+colsize);
				if (!numactive) {
					continue;
				}
				totalvar += iblsize*colsize;
			}
		}
    }
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Number of equations " << totaleq << " number of nonzero s " << totalvar;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
	
    TPZManVector<int64_t,400> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<TVar> EqValue(totalvar);
    
    //lambda for avoid repeating code
    auto AddColEqs =
        [this,&EqCol,&EqValue,&pos](const int colsize, const int colpos, const int64_t minrow)
        {
            TPZManVector<int64_t> destindices(colsize);
            for (int64_t i=0; i<colsize; i++) {
                destindices[i] = colpos+i;
            }
            this->fEquationFilter.Filter(destindices);
            for(auto jbleq=0; jbleq<destindices.size(); jbleq++) {
                if (destindices[jbleq] >= minrow) {  // Only upper triangular
                    EqCol[pos] = destindices[jbleq];
                    EqValue[pos] = 0.;
                    pos++;
                }
            }
        };

    for(auto i=0;i<nblock;i++){
		const int64_t iblsize = this->fMesh->Block().Size(i);
		const int64_t iblpos = this->fMesh->Block().Position(i);
        const int64_t numactive =
            this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t ij=0; ij<iblsize; ij++) {
            rowdestindices[ij] = iblpos+ij;
        }
        this->fEquationFilter.Filter(rowdestindices);
        // working equation by equation
        // rowdestindices contains the equation number of each element in the block number "i"
		for(auto ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            if (rowdestindices[ibleq] != ieq) {
                DebugStop();
            }
			Eq[ieq] = pos;
            //we need to add the diagonal block (only upper triangular part)
            AddColEqs(this->fMesh->Block().Size(i),
                      this->fMesh->Block().Position(i), ieq);
			
			const int64_t icfirst = nodegraphindex[i];
			const int64_t iclast = nodegraphindex[i+1];
			for(auto j=icfirst;j<iclast;j++)
            {
                // col is the block linked to block "i"
				const int64_t col = nodegraph[j];
				if (col > i) {  // Only upper triangular
					const auto colsize = this->fMesh->Block().Size(col);
					const auto colpos = this->fMesh->Block().Position(col);
					const int64_t numactive =
						this->fEquationFilter.NumActive(colpos, colpos+colsize);
					if (!numactive) {
						continue;
					}
					AddColEqs(colsize, colpos, ieq);
				}
			}
            //we need to sort the entries otherwise solvers may complain
            const auto firstentry = Eq[ieq];
            const auto nentries = pos-Eq[ieq];

            std::stable_sort(&EqCol[firstentry],
                             &EqCol[firstentry]+nentries);
            
			ieq++;
		}
    }
    Eq[ieq] = pos;
    if(pos != totalvar)
    {
        PZError<<__PRETTY_FUNCTION__<<'\n';
        PZError<<"nvars: "<<totalvar<<'\n';
        PZError<<"l pos: "<<pos<<'\n';
        DebugStop();
    }
    mat->SetData(std::move(Eq),std::move(EqCol),std::move(EqValue));
    
    return mat;
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"
#include "TPZStructMatrixOMPorTBB.h"

template class TPZSSpStructMatrixSpooles<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSSpStructMatrixSpooles<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSSpStructMatrixSpooles<STATE,TPZStructMatrixTBBFlow<STATE>>;
template class TPZSSpStructMatrixSpooles<STATE,TPZStructMatrixOMPorTBB<STATE>>;
