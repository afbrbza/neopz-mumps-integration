/**
 * @file
 * @brief Contains the implementation of TPZSSpStructMatrixMumps methods.
 * This is a specialized version that uses MUMPS solver.
 */

#include "TPZSSpStructMatrixMumps.h"
#include "TPZRenumbering.h"
#include "TPZSYSMPMatrix.h"
#include "TPZTimer.h"
#include "pzcmesh.h"
#include "pzlog.h"

#include "TPZSYSMPMumps.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar, class TPar>
TPZSSpStructMatrixMumps<TVar, TPar> * TPZSSpStructMatrixMumps<TVar,TPar>::Clone(){
    return new TPZSSpStructMatrixMumps(*this);
}

template <class TVar, class TPar>
TPZMatrix<TVar> *TPZSSpStructMatrixMumps<TVar, TPar>::Create() {
  int64_t neq = this->fMesh->NEquations();

  // if(this->fMesh->FatherMesh()) {
  //		cout << "TPZSSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
  //        DebugStop(); // WHY?
  // }

  /**
   *Longhin implementation
   */
  TPZStack<int64_t> elgraph;
  TPZVec<int64_t> elgraphindex;
  //    int nnodes = 0;
  this->fMesh->ComputeElGraph(elgraph, elgraphindex);

  TPZMatrix<TVar> *mat = SetupMatrixData(elgraph, elgraphindex);
  return mat;
}

/**
 * Implementação específica do MUMPS aqui
 * Você pode sobrescrever apenas o método SetupMatrixData se precisar
 * de um comportamento diferente para MUMPS
 */

template <class TVar, class TPar>
TPZMatrix<TVar> *TPZSSpStructMatrixMumps<TVar, TPar>::SetupMatrixData(
    TPZStack<int64_t> &elgraph,
    TPZVec<int64_t> &elgraphindex) {


  const int64_t neq = this->fEquationFilter.NActiveEquations();
  TPZSYsmpMatrixMumps<TVar> *mat = new TPZSYsmpMatrixMumps<TVar>(neq, neq);

  /**Creates a element graph*/
  TPZRenumbering graph;
  graph.SetElementsNodes(elgraphindex.NElements() - 1, this->fMesh->NIndependentConnects());

  TPZManVector<int64_t> nodegraph;
  TPZManVector<int64_t> nodegraphindex;
  /**
   *converts an element graph structure into a node graph structure
   *those vectors have size ZERO !!!
   */
  graph.ConvertGraph(elgraph, elgraphindex, nodegraph, nodegraphindex);
  /**vector sizes*/
  const int64_t nblock = nodegraphindex.NElements() - 1;
  // number of values in the sparse matrix
  int64_t totalvar = 0;
  // number of equations
  int64_t totaleq = 0;
  for (auto i = 0; i < nblock; i++) {
    const int64_t iblsize = this->fMesh->Block().Size(i);
    const int64_t iblpos = this->fMesh->Block().Position(i);
    const int64_t numactive = this->fEquationFilter.NumActive(iblpos, iblpos + iblsize);
    if (!numactive) {
      continue;
    }
    totaleq += iblsize;
    const int64_t icfirst = nodegraphindex[i];
    const int64_t iclast = nodegraphindex[i + 1];
    // longhin
    totalvar += (iblsize * (iblsize + 1)) / 2;
    for (auto j = icfirst; j < iclast; j++) {
      const int64_t col = nodegraph[j];
      if (col < i) {
        continue;
      }

      if (col == i) {
        DebugStop();
      }

      const int64_t colsize = this->fMesh->Block().Size(col);
      const int64_t colpos = this->fMesh->Block().Position(col);
      const int64_t numactive = this->fEquationFilter.NumActive(colpos, colpos + colsize);
      if (!numactive) {
        continue;
      }
      totalvar += iblsize * colsize;
    }
  }

  int64_t ieq = 0;
  // pos is the position where we will put the column value
  int64_t pos = 0;

  TPZVec<int64_t> Eq(totaleq + 1);
  TPZVec<int64_t> EqCol(totalvar);
  TPZVec<TVar> EqValue(totalvar, 0.);
  // lambda for avoid repeating code
  // lambda for avoid repeating code
  auto AddColEqs =
      [this, &EqCol, &EqValue, &pos](const int colsize, const int colpos, const int ieq) {
        TPZManVector<int64_t> destindices(colsize);
        for (int64_t i = 0; i < colsize; i++) {
          destindices[i] = colpos + i;
        }
        this->fEquationFilter.Filter(destindices);
        for (auto jbleq = 0; jbleq < destindices.size(); jbleq++) {
          const int64_t jeq = destindices[jbleq];
          if (jeq < ieq) {
            continue;
          }
          EqCol[pos] = destindices[jbleq];
          EqValue[pos] = 0.;
          pos++;
        }
      };

  for (auto i = 0; i < nblock; i++) {
    const int64_t iblsize = this->fMesh->Block().Size(i);
    const int64_t iblpos = this->fMesh->Block().Position(i);
    const int64_t numactive =
        this->fEquationFilter.NumActive(iblpos, iblpos + iblsize);
    if (!numactive) {
      continue;
    }
    TPZManVector<int64_t> rowdestindices(iblsize);
    for (int64_t i = 0; i < iblsize; i++) {
      rowdestindices[i] = iblpos + i;
    }
    this->fEquationFilter.Filter(rowdestindices);
    // working equation by equation
    for (auto ibleq = 0; ibleq < rowdestindices.size(); ibleq++) {
      if (rowdestindices[ibleq] != ieq) {
        DebugStop();
      }
      Eq[ieq] = pos;
      bool diagonalinsert = false;
      const int64_t icfirst = nodegraphindex[i];
      const int64_t iclast = nodegraphindex[i + 1];
      for (auto j = icfirst; j < iclast; j++) {
        const int64_t col = nodegraph[j];
        if (col < i) {
          continue;
        }
        // force the diagonal block to be inserted
        // the nodegraph does not contain the pointer to itself
        if (!diagonalinsert && col > i) {
          diagonalinsert = true;
          const auto colsize = this->fMesh->Block().Size(i);
          const auto colpos = this->fMesh->Block().Position(i);
          AddColEqs(colsize, colpos, ieq);
        }
        const auto colsize = this->fMesh->Block().Size(col);
        const auto colpos = this->fMesh->Block().Position(col);
        if (this->fEquationFilter.NumActive(colpos, colpos + colsize) == 0) {
          continue;
        }
        AddColEqs(colsize, colpos, ieq);
      }
      // all elements are below (last block certainly)
      if (!diagonalinsert) {
        diagonalinsert = true;
        const auto colsize = this->fMesh->Block().Size(i);
        const auto colpos = this->fMesh->Block().Position(i);
        AddColEqs(colsize, colpos, ieq);
      }
      ieq++;
    }
  }

  Eq[ieq] = pos;
  mat->SetData(std::move(Eq), std::move(EqCol), std::move(EqValue));
  return mat;
}

// Instantiate template classes for MUMPS versions
#include "TPZStructMatrixOMPorTBB.h"
#include "pzstrmatrixflowtbb.h"
#include "pzstrmatrixot.h"

template class TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOR<STATE>>;
template class TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOT<STATE>>;
template class TPZSSpStructMatrixMumps<STATE, TPZStructMatrixTBBFlow<STATE>>;
template class TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOMPorTBB<STATE>>;
