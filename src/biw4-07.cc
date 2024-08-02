// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
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

// For Plotting in Paraview
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>


using namespace Dune;

// Now a little generic programming:
// We want something derived from LocalView, without doing a virtual function lookup
// therefore we do a template
// LocalView: LocalView of Geometry
// 
template<class LocalView>
void assembleElementStiffnessMatrix(const LocalView& localView,
Matrix<double>& elementMatrix)
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
    for (size_t i=0; i<gradients.size(); i++)
      jacobianInverseTransposed.mv(referenceGradients[i][0], gradients[i]);

    for (size_t i=0; i<displacementLocalFiniteElement.size(); i++)
      for (size_t j=0; j<displacementLocalFiniteElement.size(); j++ )
        for (size_t k=0; k<dim; k++)
        {
          // We traverse the tree, get the $k$th-child an get their indezes
          // Question: Does this only work, because the function space basis is the same for each subspace?
          size_t row = localView.tree().child(k).localIndex(i);
          size_t col = localView.tree().child(k).localIndex(j);

          // here is the bilinear form of the differential equation!
          elementMatrix[row][col] += ( gradients[i] * gradients[j] )
          * quadPoint.weight() * integrationElement;
        }
    
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
template<class Basis, class Matrix>
void assembleStiffnessMatrix(const Basis& basis, Matrix& matrix)
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
    Dune::Matrix<double> elementMatrix;
    assembleElementStiffnessMatrix(localView, elementMatrix);
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

  // With this the Grid is Read and the geometry is setup.
  // We use YaspGrid for starters
  using Grid = Dune::YaspGrid<2>;
  Grid grid{ {1.0, 1.0}, {4, 4} };

  auto gv = grid.leafGridView();

  using GridView = typename Grid::LeafGridView;
  GridView gridView = grid.leafGridView();

  using namespace Functions::BasisFactory;

  // Linear Ansatz Order
  constexpr int p = 1;

  // constructing the basis:
  auto lagrangeBasis = makeBasis(gridView,
                                   power<dim>(
                                   lagrange<p>(),
                                   blockedInterleaved()));
  

  // Multiindex: Nodes x Num Dimension
  using DisplacementVector = BlockVector<FieldVector<double,dim>>;
  // Double the Multiindex: One Block for each Node of Dimension (dim x dim)
  using Matrix = BCRSMatrix<FieldMatrix<double,dim,dim>>;

  // Declare the type of the rhs
  DisplacementVector rhs;

  // initialize the datastructure
  auto rhsBackend = Functions::istlVectorBackend(rhs);
  rhsBackend.resize(lagrangeBasis); // and allocate memory for it.

  // Initialization of it's entries
  rhs = 0;

  // Declare the type of the stiffness matrix
  Matrix stiffnessMatrix;

  std::cout << "Assembling step." << std::endl;

  assembleStiffnessMatrix(lagrangeBasis, stiffnessMatrix);


  std::cout << "Dirichlet step." << std::endl;
  ////
  // Incorporate Dirichlet values.
  ////
  using DisplacementBitVector = std::vector<std::array<char,dim> >;
  DisplacementBitVector isBoundary;

  auto isBoundaryBackend = Functions::istlVectorBackend(isBoundary);
  isBoundaryBackend.resize(lagrangeBasis);

  using namespace Indices;
  for (auto&& b0i : isBoundary)
    for (std::size_t j=0; j<b0i.size(); ++j)
      b0i[j] = false;

  ////
  // Convenience Function for Boundary DOFs
  Functions::forEachBoundaryDOF(
    lagrangeBasis,
    [&] (auto&& index) {
      isBoundaryBackend[index] = true;
    });
  
  ////
  // Genau den hier machen
  using Coordinate = GridView::Codim<0> ::Geometry::GlobalCoordinate;
  using DisplacementRange = FieldVector<double,dim>;

  // Funktion die den Dirichlet-Rand vorgibt
  auto&& g = [](Coordinate x)
    {
      return DisplacementRange{0.0, (x[1] - 1.0 < 1e-8) ? 1.0 : 0.0};
    };
  Functions::interpolate(lagrangeBasis,
      rhs,
      g,
      isBoundary);

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

  // that way the Dirichlet entries are already correct.
  DisplacementVector x = rhs;
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

  using VelocityRange = FieldVector<double,dim>;
  
  auto displacementFunction
    = Functions::makeDiscreteGlobalBasisFunction<DisplacementRange>(lagrangeBasis, x);

  SubsamplingVTKWriter<GridView> vtkWriter(
    gridView,
    refinementLevels(2));

  vtkWriter.addVertexData(
    displacementFunction,
    VTK::FieldInfo("displacement", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.write("displacement-result");
}
