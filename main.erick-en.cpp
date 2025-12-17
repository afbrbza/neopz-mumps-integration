#include "TPZMultiphysicsCompMesh.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSimpleTimer.h"
#include "TPZStructMatrixOMPorTBB.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzfmatrix.h"
#include "pzlog.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include <DarcyFlow/TPZDarcyFlow.h>
#include <TPZGenGrid2D.h>
#include <TPZGenGrid3D.h>
#include <TPZLinearAnalysis.h>
#include <iostream>

using namespace std;
/* just a facilitation for the material ids I will use this in the boundary conditions too
It's more of a parametrization of things */
enum EnumMatids { EMatId = 1,
                  EBottom = 2,
                  ETop = 3,
                  ELeft = 4,
                  ERight = 5 };
/* We are passing to createMeshWithgenGrid a vector with the number of divisions and the minimum and maximum coordinates in x.
In the generator, we set the element type */
// ... (suas funções de mesh) ...
TPZGeoMesh *createMeshWithGenGrid(const TPZVec<int> &nelDiv, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX) {
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  // creating a 2d object named generator that will use a number of divisions (nelDiv) in x and y from min to max
  TPZGenGrid2D generator(nelDiv, minX, maxX);
  // creating the quadrilateral geometric mesh. Must specify the variable type (MMeshType)
  generator.SetElementType(MMeshType::EQuadrilateral);
  /* here it creates the mesh by passing to gmesh the configurations we passed and the boundary conditions
  First, we tell it to read the material ID. Then we give */
  generator.Read(gmesh, EMatId);
  generator.SetBC(gmesh, 4, EBottom); // bottom
  generator.SetBC(gmesh, 5, ERight);  // right
  generator.SetBC(gmesh, 6, ETop);    // top
  generator.SetBC(gmesh, 7, ELeft);   // left
  return gmesh;
}
TPZGeoMesh *createRectangularGmesh() {
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  /* setting the coordinate vectors of the nodes because TPZ wants them to be of type double. This is not the most efficient way to do it, but the important thing is that it works for now; later, it needs to be made more presentable */

  TPZVec<REAL> coord0 = {0., 0., 0.};
  TPZVec<REAL> coord1 = {1., 0., 0.};
  TPZVec<REAL> coord2 = {1., 1., 0.};
  TPZVec<REAL> coord3 = {0., 1., 0.};
  TPZVec<REAL> coord4 = {2., 0., 0.};
  TPZVec<REAL> coord5 = {2., 1., 0.};
  /* creating the nodes with reference to the created vectors. The first entry we have is the node id, the second is the coordinates, and the third is the reference of the mesh to which we are adding the node */
  TPZGeoNode nod0(0, coord0, *gmesh);
  TPZGeoNode nod1(1, coord1, *gmesh);
  TPZGeoNode nod2(2, coord2, *gmesh);
  TPZGeoNode nod3(3, coord3, *gmesh);
  TPZGeoNode nod4(4, coord4, *gmesh);
  TPZGeoNode nod5(5, coord5, *gmesh);
  /* then we add each of the nodes created in gmesh because apparently they are not added automatically
  after its creation*/
  /* First, it gets the node vector of my geometric mesh (pointer to the obj) passes it to NodeVec which returns my node vector, and we resize it to a vector of 4 positions (# I don't know why). Then it adds the respective nodes to positions 0, 1, 2, 3. Note: in PZ, index is the position in the element vector, and Id is any name we try to make identifiable with the index */
  gmesh->NodeVec().Resize(6);
  gmesh->NodeVec()[0] = nod0;
  gmesh->NodeVec()[1] = nod1;
  gmesh->NodeVec()[2] = nod2;
  gmesh->NodeVec()[3] = nod3;
  gmesh->NodeVec()[4] = nod4;
  gmesh->NodeVec()[5] = nod5;
  /* So here we create a vector of Node Indexes which will be used to create the quadrilateral element (Ps: remember to identify index with id so as not to get lost and be able to identify each node in the element more easily)*/
  //====================== First Element======================
  TPZManVector<int64_t, 4> nodeIndexes(4);
  nodeIndexes[0] = 0;
  nodeIndexes[1] = 1;
  nodeIndexes[2] = 2;
  nodeIndexes[3] = 3;
  int matid = 1;     // pass the ids as integer identifiers
  int64_t index = 0; // pass the ids as integer64 identifiers
  TPZGeoEl *gel = gmesh->CreateGeoElement(EQuadrilateral, nodeIndexes, EMatId, index);
  /* So we create a geometric element of type Equadrilateral, which is an Enum. A geometric element is an element that contains only information about the mesh's geometry, such as size, shape, divisions, etc., which are defined in PZ with the node indexes we created right above */
  /* The matId is the equation, it is the weak form of the differential equation
  related to this element. This has ID 1, so when I create the material, it must have the same ID (1), otherwise, things go wrong, or maybe not. For a heterogeneous element, we can set different materials for each part of the mesh, whether defined in the same computational mesh or not. (I'll explain what that is shortly...) Note: each element is associated with a weak formulation, which is a reconfiguration of the original differential equation with mathematical tricks - like substitutions, derivatives, and integrals to make it easier to
  solve. As soon as you solve this little equation, you can find a solution for the original monster. That's more or less it. You should look up the Galerkin method, it's about this. */
  //====================== Second Element======================
  nodeIndexes[0] = 1;
  nodeIndexes[1] = 4;
  nodeIndexes[2] = 5;
  nodeIndexes[3] = 2;
  index = 1;
  TPZGeoEl *gel2 = gmesh->CreateGeoElement(EQuadrilateral, nodeIndexes, EMatId, index);
  /* ====================== Creating Boundary Conditions ======================
  Basically what happens to the element in some places where we know. For example
  in an embedded bar */
  // int idriri = 2;
  gel->CreateBCGeoEl(4, EBottom);
  gel2->CreateBCGeoEl(4, EBottom);
  gel2->CreateBCGeoEl(5, ERight);
  gel2->CreateBCGeoEl(6, ETop);
  gel->CreateBCGeoEl(6, ETop);
  gel->CreateBCGeoEl(7, ELeft);
  /* this element is of the type:
  _____6_________6________
  | | |
  | (0) 5| (1) | 5
  7 | |7 |
  |__________|___________|
  4 4
  */
  // gmesh -> Print(std::cout);
  /* for it to print the mesh before, you need to set the neighbors so that
  the mesh conforms. This way, you connect the things - points and elements - so that this tangle of things makes sense as one thing (for example, so that bars connect forming a coherent structure)*/
  gmesh->SetDimension(2);
  gmesh->BuildConnectivity();
  /* Note: until then this element is not GeoBlent, meaning it does not have curved mapping */
  /* This command creates the boundary condition elements (BC) associated with the geometric element. It requires some information such as side and BCs
  both as integers (int).
  This method creates a geometric element adjacent to another element. This must maintain an order. Here I am setting the sides of the two elements counter-clockwise, so side 4 is South, side 5 is East, side 6 is North, and side 7 is West (as seen above). And as you can see, the East sides (5) of element 0 and West sides (7) of element 1 are joined on the inside, so BCs are not set for them*/
  return gmesh;
}
/* All of this created so far is the geometric mesh without any attributes. */
/*
======================= COMPUTATIONAL MESH ============================
From here we will create a computational mesh. What is that?
This is the mesh set with the physical and material attributes of the gmesh.
For this, you first have to create a geometric mesh - you will see that soon all cmeshes are created from a gmesh. Another important thing is that if your material is heterogeneous, that is, has different parts with different materials, you can set this in the computational mesh, or depending on the case
a computational mesh with different materials for each part of the geometric mesh.
As for the weak formulation, this will depend on the material chosen for the mesh, which I already explained how it works more or less above, just take a look (CTRL F "weak")*/
TPZCompMesh *createCompMesh(TPZGeoMesh *gmesh) {
  /* first we create a CompMesh, which I don't know what exactly it is, and we set a pointer to return it and base this new element on a gmesh*/
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  /* We set the dimension of the model as equal to the dimension of our gmesh model, which for now is 2 */
  cmesh->SetDimModel(gmesh->Dimension());
  // Now we set the polynomial order of approximation
  cmesh->SetDefaultOrder(1);
  /* this is basically the polynomial order of the shape functions used in the approximation, which are formed with hat functions */
  /* And now we set the type of approximation space we want to use, for example H1 Hdiv etc. I don't know the difference yet, but we will use H1 here
  From what I understood, these would be traditional FEM with hat functions that are kind of functions used as a basis that may or may not be linear (in general they are) that have their sample space (let's say) varying between 0 and 1, meaning being a function of order 1 it does this:
  (y)
  y=1 ____^
  | /\ /\
  | / \ / \
  | / \ / \ ...
  |---------------------> (x) */
  cmesh->SetAllCreateFunctionsContinuous(); // H1
  /*====================== Creating the Material and the Weak Formulation ======================
  then we will add the materials and the weak formulations
  In PZ, from what I understood, these material functions already come with the weak formulation of the thing, so it would only be a matter of calling and setting them correctly*/
  int matId = 1; // must be the same id as the geometric element
  TPZDarcyFlow *mat = new TPZDarcyFlow(EMatId, gmesh->Dimension());
  /* This basically defines that our problem is a Darcy flow, so it goes to the Darcy library and uses its formulations to solve the problem. CTRL it if you want to see how it works*/
  /* we set the permeability value equal to 1.0, that is, the same for the entire domain
  and equal to an ideal fluid*/
  mat->SetConstantPermeability(1.0);
  // we insert this material into the computational mesh
  cmesh->InsertMaterialObject(mat);
  /* So far, practically everything has been set, but nothing has been built yet. Recapping, we are setting here that for every geometric element with Id 1 that has matId being 1, it says that the thing will have the same material as the one defined above. Is that clear? I don't know. More or less*/
  /*============boundary conditions====================*/
  int dirBC = 2;
  /* To create boundary conditions, the function in PZ is based on val1 and val2, where val1 is a matrix that goes to the stiffness matrix relative to the boundary condition, and val2 is the load vector of the boundary condition [K].{x} = {F}. What is a Neumann condition in the traditional sense? I have no idea, but when we have this, everything goes to the load vector {F}, so val1 = 0 and val2 is what matters. The two are only useful in a mixed condition (which I also don't know what it is). From what I understood, Neumann and Dirichlet are used for different situations, one where we have
  the boundary conditions and the other when we have some variable. For example, in an embedded bar, we can have the boundary condition that the displacement is zero at the constraint and work with that, or even say that there is a load on the bar and work with that. In a mixed regime, we can work with both things.*/
  TPZManVector<REAL, 1> val2(1, 3.); // Parts that goes to the RHS
  TPZFMatrix<REAL> val1(1, 1, 0.);   // Parts that goes to the Stiffness matrix
  int diritype = 0, neumanntype = 1, mixedtype = 2;
  /* then we create the two boundary conditions that we will use.
  Note: the REAL in the Object is because of Fran who did these things with complex numbers
  complex (an electric folks thing)
  The createBC works as follows: here we create the condition name by passing the material to the function, in it we set the material (a bit redundant but easily modifiable) we pass the ID of the boundary condition, right after we have the type, we have 3 as you can see above*/
  /* we insert the 4 boundary conditions into the computational mesh.
  The inputs for createBC are material, where the boundary condition is, its initial and final values*/
  val2[0] = 0.;
  TPZBndCondT<REAL> *bcond = mat->CreateBC(mat, EBottom, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  val2[0] = 1.;
  bcond = mat->CreateBC(mat, ERight, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  val2[0] = 0.;
  bcond = mat->CreateBC(mat, ETop, neumanntype, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  val2[0] = 3.;
  bcond = mat->CreateBC(mat, ELeft, diritype, val1, val2);
  cmesh->InsertMaterialObject(bcond);
  /* what remains now is to do the autobuild which will build everything for us and return our mesh*/
  cmesh->AutoBuild();
  return cmesh;
}
// ... (suas funções de mesh) ...
int main(int argc, char const *argv[]) {
#ifdef PZ_LOG
  // Initialize PZLOG
  TPZLogger::InitializePZLOG();
#endif
  /* number of threads for the skyline matrix */
  const int nthreads = 16;
  /* The number of threads basically tells how you want this problem to run. It refers to the processor speed. Two means it will run with two threads, i.e., two tasks simultaneously, and so on. */
  // ====================== Creating the Geometric Mesh ======================
  const bool isUseGenGrid = true;
  TPZGeoMesh *gmesh = nullptr;
  if (isUseGenGrid) {
    const int neldiv = 2;
    gmesh = createMeshWithGenGrid({neldiv, neldiv}, {0., 0.}, {2., 1.});
  } else {
    gmesh = createRectangularGmesh();
  }
  gmesh->Print(std::cout);
  /* To see this thing in Paraview, it seems like you just have to call this function:*/
  std::ofstream out("gmesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  // ====================== Creating the Computational Mesh ======================
  TPZCompMesh *cmesh = createCompMesh(gmesh);
  cmesh->Print(std::cout);
  // ====================== Linear Analysis ======================
  /* Now we are going to do the local matrices and the assembly into the global matrix. Solve the system and plot this thing in vtk paraview. There is a class that already automates this for us, TPZAnalyses. We will use its daughter for linear analysis */
  /* it receives the computational mesh in the constructor. The second parameter is for whether or not to reorder the equations (I think). Let's set it to true for now*/
  TPZLinearAnalysis an(cmesh);
// #if ${USE_MKL ^ ^} == USE_MKL
#if USE_MKL == USE_MKL
  // MKL version
  /* So we have to choose which type of matrix to use. Here I'm going to try to do it with the paradiso matrix which is the best (which apple cannot use) */
  TPZSSpStructMatrix > matsp(cmesh);
  // setting the number of threads for the matrix
  matsp.SetNumThreads(nthreads);
  /* telling the analysis to use this */
  an.SetStructuralMatrix(matsp);
  // process to determine the solver to be used with this type of matrix
  TPZStepSolver step;
  // step.SetDirect(ECholesky); // commented because it does not work with MKL, the solver itself defines it based on the matrix type
  an.SetSolver(step);
  an.Assemble();
  auto mCast = an.MatrixSolver().Matrix();
  mCast->SetDefPositive(true);
  an.Solve();
#else
  // non-MKL version
  /* Just to follow the tutorial, we'll do it the way Nathan did with a worse solver. The skyline counts the last non-null position for each column and stores everything somewhere I didn't understand (have to check this later), but according to Nathan, it looks like a building profile (I don't know what that means)*/
  TPZSkylineStructMatrix matsp(cmesh);
  // setting the number of threads for the matrix
  matsp.SetNumThreads(nthreads);
  /* telling the analysis to use this */
  an.SetStructuralMatrix(matsp);
  // process to determine the solver to be used with this type of matrix
  TPZStepSolver step;
  /* using a direct solver which is the opposite of the iterative mode (I also don't know what that means) and uses ECholesky (which I also don't know what it is, but I believe it is a definition of positive definite and symmetric)*/
  step.SetDirect(ECholesky);
  // using this solver for the analysis
  an.SetSolver(step);
  an.Run(); // this just runs the whole analysis
#endif
  an.Solution().Print("Solution = "); /// just to show the solution in the terminal
  //====================== Post Processing ======================
  const std::string plotfile = "postproct";
  /* A string with the name of the thing must be without the ".vtk" because it adds it later*/
  constexpr int vtkRes{0};
  /* from what I understood, vtk has some nodal solutions, and if the element was quadratic, it kind of loses the connections. So Phil refines the elements and creates a node between edges and samples there too, which is basically adding intermediate nodes
  to maintain the connection.
  Nathan says it has quadratic support, but only Lagrangian squares. What is that?
  _________________________________________________________________________________________
  A Lagrangian square is a quadratic form associated with the method of Lagrange multipliers. In simple terms, it is the mathematical expression obtained by applying the Lagrange method to solve optimization problems with constraints, where the objective function and the constraints are combined into a single function called the Lagrangian.
  🔎 Context
  In optimization problems, we want to maximize or minimize a function f(x) subject to certain constraints g_i(x)=0.
  The Lagrange method constructs the function:
  The so-called Lagrangian square appears when analyzing the quadratic form associated with the second derivative (Hessian) of the Lagrangian, used to verify whether a critical point is a minimum, maximum, or saddle point.
  📌 In practical terms
  The Lagrangian square is a quadratic expression that evaluates the curvature of the objective function taking into account the constraints:
  It is used to determine the nature of the solution found:
  If the Lagrangian square is positive definite, the point is a local minimum.
  If it is negative definite, the point is a local maximum.
  If it is indefinite, the point is a saddle point.
  In summary: the Lagrangian square is the quadratic form derived from the Lagrangian, used to analyze the nature of critical points in constrained optimization problems. It is an essential tool to confirm whether the obtained solution is a minimum, maximum, or saddle point.*/
  /* So we create a vector of strings with what we want to post-process */
  TPZManVector fields = {"Pressure", "Flux"};
  /* Next, we create an auto object with TPZVTKGenerator which receives the meshes, the fields I want to post-process, and the resolution*/
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.Do();
  return 0;
}
/*===============================================================================================
WHAT IS A CONNECT?
A Connect is a data structure that represents a set of shape functions.
This is what PZ printed from the computational mesh. Let's see what a connect is:
Index 0 TPZConnect : Sequence number = 0 Order = 1 NState = 1 NShape 1 IsCondensed 0 IsLagrMult 0
Equation = 0 NumElCon = 3 Block size 1 Solution 0
# Sequence number = 0 -> where it goes in the global matrix, its order in the matrix
Correction: it doesn't go into the global matrix, it is a unique identifier of the connect
# Order = 1 -> polynomial order of the connect, meaning you can have connects of order 1, 2, 3 etc
in selected connects and not in the entire mesh
# NState = 1 -> number of state variables, meaning it stores how many variables
different you want to associate with this connect, which are the physical behavior of the thing in a given direction.
# NShape 1 -> number of shape functions associated with this connect
# IsCondensed 0 -> if it is related to static condensation (I don't know what that is)
# IsLagrMult 0 -> if it is related to Lagrange multipliers (I don't know what that is, but they say it's more for mixed formulations)
# Equation = 0 -> I think it's which equations. I'm not sure Nathan doesn't remember
But according to Phil, this is the number that associates the global matrix
# NumElCon = 3 -> number of elements associated with this connect, including itself
# Block size 1 -> block size, I'm not sure what it is
# Solution 0 -> multiplier of the shape function associated with this connect
================================================================================================
TYPES OF MATRICES USED IN PZ
CHOLESKY - used for symmetric and positive definite systems.
"A symmetric matrix is one that satisfies the condition A = A^T. while a positive definite matrix
positive definite is one that satisfies V (this means for every single one, the symbol is actually an upside-down A) x (vector) =! 0, x^t @ A @ x > 0
these matrices can be written as:
A = L @ L^T "
SKYLINE - used for symmetric and positive definite systems
SSPARSE - used for sparse non-symmetric systems
SSSPARSE - used for sparse symmetric systems
"A sparse matrix is one that has a large quantity of elements equal to zero."
LDLT - used for symmetric and indefinite systems
A = L @ D @ L^t
ELU - used for non-symmetric systems
===============================================================================================
*/