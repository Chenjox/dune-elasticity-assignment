// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// this is first
#include "config.h"
#include "dune/biw4-07/biw4-07.hh"
#include "dune/common/fmatrix.hh"
#include <cstddef>
#include <vector>
//#ifdef HAVE_CONFIG_H
//#endif
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/yaspgrid.hh>
// With this the Mesh will be isoparametric for 
//#include <dune/foamgrid/foamgrid.hh> // <-- does only work for quads
//#include <dune/grid/uggrid.hh>
//#include <dune/curvedgrid/curvedgrid.hh>
//#include "dune/gmsh4/gmsh4reader.hh"
//#include <dune/gmsh4/gmsh4reader.impl.hh>

// Matrix types and storage, as well as convenience
#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// Now we will use Dune Functions
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

// These are functions defined on a grid
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>

// For Plotting in Paraview
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>


using namespace Dune;

// Now a little generic programming:
// We want something derived from LocalView, without doing a virtual function lookup
// therefore we do a template
// LocalView: LocalView of Geometry
// 
template<class LocalView, class GridFunction>
void assembleElementStiffnessMatrix(
  const LocalView& localView,
  Matrix<double>& elementMatrix,
  const GridFunction& displacements)
{
  // what type of element do I have?
  using Element = typename LocalView::Element;
  const Element element = localView.element();
  constexpr int dim = Element::dimension;

  // take the geometry information. TODO: CurvedGrid for Isoparametric Elements.
  auto geometry = element.geometry();

  // Set the size
  elementMatrix.setSize(localView.size(), localView.size());
  // Set the numbers
  elementMatrix = 0;


  using namespace Indices;
  const auto& displacementLocalFiniteElement
    = localView.tree().child(0).finiteElement();

  // How many Dofs, on the finite element in total?
  int num_dofs = displacementLocalFiniteElement.localBasis().size();

  // take necessary integration order.
  int order = 2*(dim*displacementLocalFiniteElement.localBasis().order());
  const auto& quad = QuadratureRules<double,dim>::rule(element.type(), order);

  // Loop over all quadrature points
  for (const auto& quadPoint : quad)
  {
    // The transposed inverse Jacobian of the map from the
    // reference element to the element
    const auto jacobianInverseTransposed
    = geometry.jacobianInverseTransposed(quadPoint.position());
    // The multiplicative factor in the integral transformation formula (chain-rule)
    const auto integrationElement
    = geometry.integrationElement(quadPoint.position());

    // 
    std::vector<FieldMatrix<double,1,dim> > referenceGradients;
    displacementLocalFiniteElement.localBasis().evaluateJacobian(quadPoint.position(),referenceGradients);
    // Compute the shape function gradients on the grid element
    std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
    for (size_t i=0; i<gradients.size(); i++){
      gradients[i] = 0.0;
      jacobianInverseTransposed.mv(referenceGradients[i][0], gradients[i]);
    }

    // A two dimensional list of FieldMatrizes, 
    // storing the strains of every shape function of every strain.
    // seems inefficient to me, but it works.
    std::vector<std::array<FieldMatrix<double, dim, dim>,dim>> deltaLinStrain(num_dofs);
    // Loop over the Dofs
    for (size_t i=0; i<num_dofs; i++)
    {
      for (size_t k=0; k<dim; k++)
      {
        // 
        FieldMatrix<double,dim,dim> deformationGradient(0);
        deformationGradient[k] = gradients[i];
        //elementMatrix[row][col] += ( gradients[i] * gradients[j] )* quadPoint.weight() * integrationElement;

        FieldMatrix<double, dim, dim> linearisedStrains(0);
        for (size_t m=0; m<dim; m++)
          for (size_t n=0; n<dim; n++){
            linearisedStrains[m][n] += 0.5 * (deformationGradient[m][n] + deformationGradient[m][n]);
        }
      }
    }

    // get the _real_ deformation gradient
    typename GridFunction::DerivativeType localConfGrad;
    if (displacements.isDefinedOn(element))
        displacements.evaluateDerivativeLocal(element, quadPoint, localConfGrad);
    else
      displacements.evaluateDerivative(geometry.global(quadPoint),localConfGrad);

  }
}

template<class Basis, class Matrix>
void setOccupationPattern(const Basis& basis, Matrix& matrix)
{
  enum {dim = Basis::GridView::dimension};
  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet nb;
  nb.resize(basis.size(), basis.size());
  // A view on the FE basis on a single element
  auto localView = basis.localView();
  // Loop over all leaf elements
  for(const auto& element : elements(basis.gridView()))
  {
    // Bind the local view to the current element
    localView.bind(element);
    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<localView.size(); i++)
    {
      // Global index of the i-th local degree of freedom of the current element
      auto row = localView.index(i);
      for (size_t j=0; j<localView.size(); j++ )
      {
          // Global index of the j-th local degree of freedom of the current element
        auto col = localView.index(j);
        nb.add(row[0],col[0]); // im ersten Multiindex steht der block index!
      }
    }
  }
  // Give the matrix the occupation pattern we want.
  using namespace Indices;
  nb.exportIdx(matrix);
}

template<class Matrix, class MultiIndex>
decltype(auto) matrixEntry(
Matrix& matrix, const MultiIndex& row, const MultiIndex& col)
{
  using namespace Indices;
  return matrix[row[0]][col[0]][row[1]][col[1]];
}

// Assemble the Laplace stiffness matrix on the given grid view
template<class Basis, class Matrix, class GridFunction>
void assembleStiffnessMatrix(const Basis& basis, Matrix& matrix, GridFunction& displacements)
{
  // Set matrix size and occupation pattern
  setOccupationPattern(basis, matrix);
  // Set all entries to zero
  matrix = 0;
  // A view on the FE basis on a single element
  auto localView = basis.localView();
  // A loop over all elements of the grid
  for (const auto& element : elements(basis.gridView()))
  {
    // Bind the local FE basis view to the current element
    localView.bind(element);

    // Now let’s get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    // displacement function must derivative from
    auto d = derivative(displacements);
    Dune::Matrix<double> elementMatrix;
    assembleElementStiffnessMatrix(localView, elementMatrix, displacements);
    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i=0; i<elementMatrix.N(); i++)
    {
    // The global index of the i-th local degree of freedom
    // of the current element
    auto row = localView.index(i);
    for (size_t j=0; j<elementMatrix.M(); j++ )
      {
      // The global index of the j-th local degree of freedom
      // of the current element
        auto col = localView.index(j);
        matrixEntry(matrix, row, col) += elementMatrix[i][j];
      }
    }
  }
}


int main(int argc, char** argv)
{
  // Maybe initialize MPI
  MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  constexpr int dim = 2;


  /////////////////////////////////////////////////////////
  // GRID DATA
  /////////////////////////////////////////////////////////
  // With this the Grid is Read and the geometry is setup.
  // We use YaspGrid for starters
  using Grid = Dune::YaspGrid<2>;
  Grid grid{ {1.0, 1.0}, {4, 4} };

  //auto gv = grid.leafGridView();

  using GridView = typename Grid::LeafGridView;
  GridView gridView = grid.leafGridView();

  using namespace Functions::BasisFactory;

  using Coordinate = GridView::Codim<0> ::Geometry::GlobalCoordinate;

  //////////////////////////////////////////////////////////
  // Function Space Basis
  //////////////////////////////////////////////////////////

  // Linear Ansatz Order
  constexpr int p = 1;

  // constructing the basis:
  auto lagrangeBasis = makeBasis(gridView,
                                   power<dim>(
                                   lagrange<p>(),
                                   blockedInterleaved())); 
  // Blocked interleaved allows for efficent sparsity structure.
  

  ///////////////////////////////////////////////////////////
  // Setting up Storage and Matrizes
  ///////////////////////////////////////////////////////////

  // Every Node has `dim` Degrees of Freedom.
  using DisplacementRange = FieldVector<double,dim>;
  // Nodal Matrix
  using DisplacementMatrix = FieldMatrix<double,dim,dim>;
  // Multiindex: Nodes x Num Dimension
  using DisplacementVector = BlockVector<DisplacementRange>;
  // Double the Multiindex: One Block for each Node of Dimension (dim x dim)
  // 
  using Matrix = BCRSMatrix<DisplacementMatrix>;
  // Marking displacement
  using DisplacementBitVector = std::vector<std::array<char,dim> >;

  // Declare the Boundary Vector, and the Right-Hand Side
  DisplacementBitVector isBoundary;
  DisplacementVector rhs;
  Matrix stiffnessMatrix;

  // 'rhsBackend' now connects rhs with `dune-functions`
  auto rhsBackend = Functions::istlVectorBackend(rhs);
  rhsBackend.resize(lagrangeBasis);

  rhs = 0;
  DisplacementVector x = rhs;

  auto isBoundaryBackend = Functions::istlVectorBackend(isBoundary);
  isBoundaryBackend.resize(lagrangeBasis);

  auto displacementFunction = Functions::makeDiscreteGlobalBasisFunction<DisplacementRange>(lagrangeBasis, x);

  using namespace Indices;
  for (auto&& b0i : isBoundary)
    for (std::size_t j=0; j<b0i.size(); ++j)
      b0i[j] = false;

  // Declare the type of the stiffness matrix

  std::cout << "Assembling step." << std::endl;

  assembleStiffnessMatrix(lagrangeBasis, stiffnessMatrix, displacementFunction);

  std::cout << "Dirichlet step." << std::endl;

  ////
  // Convenience Function for Boundary DOFs
  Functions::forEachBoundaryDOF(
    lagrangeBasis,
    [&] (auto&& index) {
      isBoundaryBackend[index] = true;
    });
  
  ////
  // Genau den hier machen

  // Funktion die den Dirichlet-Rand vorgibt
  auto&& g = [](Coordinate x)
    {
      return DisplacementRange{0.0, (x[1] - 1.0 < 1e-8) ? 1.0 : 0.0};
    };
  
  // modify rhs to incorporate dirichlet values
  Functions::interpolate(lagrangeBasis,
      rhs,
      g,
      isBoundary);
      
  //stiffnessMatrix.mmv(x,rhs);
  
  auto localView = lagrangeBasis.localView();
  for(const auto& element : elements(gridView))
  {
    localView.bind(element);
    for (size_t i=0; i<localView.size(); ++i)
    {
      auto row = localView.index(i);
      // If row corresponds to a boundary entry,
      // modify it to be an identity matrix row.
      if (isBoundaryBackend[row])
        for (size_t j=0; j<localView.size(); ++j)
        {
        auto col = localView.index(j);
        matrixEntry(stiffnessMatrix, row, col) = (i==j) ? 1 : 0;
        }
    }
  }

  //////////////////////////////////////////////////
  // Solving the ~linear~ system
  //////////////////////////////////////////////////

  // Turn the matrix into a linear operator
  std::cout << "Solving step." << std::endl;

  MatrixAdapter<Matrix,DisplacementVector,DisplacementVector> stiffnessOperator(stiffnessMatrix);

  // Sequential incomplete LU decomposition as the preconditioner
  SeqILU<Matrix,DisplacementVector,DisplacementVector> preconditioner(stiffnessMatrix,1.0); // Relaxation factor
  // Preconditioned conjugate gradient solver
  CGSolver<DisplacementVector> cg(
    stiffnessOperator,
    preconditioner,
    1e-8, // Desired residual reduction factor
    50, // Maximum number of iterations
    2); // Verbosity of the solver
  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;
  // Solve!
  cg.apply(x, rhs, statistics);
  

  SubsamplingVTKWriter<GridView> vtkWriter(
    gridView,
    refinementLevels(2));

  vtkWriter.addVertexData(
    displacementFunction,
    VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.write("displacement-result");
}
