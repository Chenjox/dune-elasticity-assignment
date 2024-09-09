// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// this is first
#include "config.h"
#include "dune/common/fmatrix.hh"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>
//#ifdef HAVE_CONFIG_H
//#endif
#include "dune/biw4-07/biw4-07.hh"

#include <dune/common/exceptions.hh>         // We use exceptions
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI

//#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

// Matrix types and storage, as well as convenience
#include <dune/istl/matrix.hh>
//#include <dune/istl/vector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

// Now we will use Dune Functions
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

// These are functions defined on a grid
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/curvedgrid/gridfunctions/discretegridviewfunction.hh>
//#include "dune/gmsh4/gmsh4reader.hh"
//#include <dune/gmsh4/gmsh4reader.impl.hh>

// For Plotting in Paraview
//#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>
#include <dune/vtk/vtkwriter.hh>

using namespace Dune;

// Now a little generic programming:
// We want something derived from LocalView, without doing a virtual function
// lookup therefore we do a template LocalView: LocalView of Geometry
//
template <class LocalView, class LocalFunctionDisp, class Material>
void assembleElementStiffnessMatrix(const LocalView &localBasisView,
                                    const LocalFunctionDisp &localDisplacements,
                                    const Material &material,
                                    Matrix<double> &elementMatrix,
                                    std::vector<double> &residualVector) {
  // what type of element do I have?
  using Element = typename LocalView::Element;
  const Element element = localBasisView.element();
  constexpr int dim = Element::dimension;

  // take the geometry information, constructed from the local displacements
  auto geometry = element.geometry();

  // Set the size of the element stiffness and residual vector
  elementMatrix.setSize(localBasisView.size(), localBasisView.size());
  residualVector.resize(localBasisView.size());

  elementMatrix = 0;
  residualVector.assign(localBasisView.size(), 0.0);

  auto const &localFE = localBasisView.tree().child(0).finiteElement();
  auto const &localBasis = localFE.localBasis();

  // How many NODES, on the finite element in total?
  int num_nodes = localBasis.size();

  // take necessary integration order.
  int order = (2 * dim * localBasis.order());
  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);

  auto localDisplacmentDerivative = derivative(localDisplacements);

  // Loop over all quadrature points
  for (const auto &quadPoint : quad) {

    // std::cout << "T" << std::endl;

    // The transposed inverse Jacobian of the map from the
    // reference element to the element to the spatial coordinates
    // This transforms the scalar gradients of the shape functions, not the
    // vector gradient
    const auto jacobianInverseTransposed =
        geometry.jacobianInverseTransposed(quadPoint.position());

    // The multiplicative factor in the integral transformation formula
    // (chain-rule)
    const auto integrationElement =
        geometry.integrationElement(quadPoint.position());

    //
    std::vector<FieldMatrix<double, 1, dim>> referenceGradients;
    localBasis.evaluateJacobian(quadPoint.position(), referenceGradients);
    // Compute the shape function gradients on the grid element
    std::vector<FieldVector<double, dim>> gradients(referenceGradients.size());
    for (size_t i = 0; i < gradients.size(); i++) {
      gradients[i] = 0.0;
      jacobianInverseTransposed.mv(referenceGradients[i][0], gradients[i]);
      // std::cout << gradients[i] << std::endl;
    }

    // A two dimensional list of FieldMatrizes,
    // storing the strains of every shape function of every strain.
    // seems inefficient to me, but it works.
    std::vector<std::array<FieldMatrix<double, dim, dim>, dim>> deltaLinStrain(
        num_nodes);
    std::vector<std::array<FieldMatrix<double, dim, dim>, dim>> sortedGradients(
        num_nodes);
    // Loop over the Dofs
    for (size_t i = 0; i < num_nodes; i++) {
      for (size_t k = 0; k < dim; k++) {
        //
        FieldMatrix<double, dim, dim> displacementGradiente(0);
        displacementGradiente[k] = gradients[i];
        // elementMatrix[row][col] += ( gradients[i] * gradients[j] )*
        // quadPoint.weight() * integrationElement;

        FieldMatrix<double, dim, dim> linearisedStrains(0);

        Dune::BIW407::linearisedStrain(displacementGradiente,
                                       linearisedStrains);
        // std::cout << linearisedStrains << std::endl;

        deltaLinStrain[i][k] = linearisedStrains;
        sortedGradients[i][k] = displacementGradiente;
      }
    }

    // now to the displacement gradient in spatial coordinates
    auto displacementGradient =
        localDisplacmentDerivative(quadPoint.position());

    Dune::FieldMatrix<double, dim, dim> InverseDeformationGradient =
        displacementGradient;

    InverseDeformationGradient *= -1.0;

    InverseDeformationGradient[0][0] += 1.0;
    InverseDeformationGradient[1][1] += 1.0;

    auto deformationGradient =
        Dune::BIW407::inverse(InverseDeformationGradient);

    FieldMatrix<double, dim, dim> cauchyStresses(0);
    FieldMatrix<double, dim, dim> cauchyStressInkrement(0);

    material.cauchyStresses(deformationGradient, cauchyStresses);
    // std::cout << deformationGradient << std::endl;
    // std::cout << cauchyStresses << std::endl;

    // Geometrical Tangent!
    //
    for (int row = 0; row < num_nodes; row++) {
      for (int col = 0; col < num_nodes; col++) {
        FieldVector<double, dim> firstEnergy(0.0);
        cauchyStresses.mv(gradients[row], firstEnergy);
        auto other = gradients[col];
        double value = other.dot(firstEnergy);
        elementMatrix[dim * row][dim * col] +=
            value * quadPoint.weight() * integrationElement;
        elementMatrix[dim * row + 1][dim * col + 1] +=
            value * quadPoint.weight() * integrationElement;
      }
    }
    // Material Tangent
    for (int row = 0; row < num_nodes; row++) {
      for (int i = 0; i < dim; i++) {
        auto realDeltStrain = deltaLinStrain[row][i];
        cauchyStressInkrement = 0;
        material.cauchyStressInkrement(deformationGradient, realDeltStrain,
                                       cauchyStressInkrement);
        for (int col = 0; col < num_nodes; col++) {
          for (int j = 0; j < dim; j++) {
            auto virtDeltStrain = deltaLinStrain[col][j];

            elementMatrix[dim * row + i][dim * col + j] +=
                Dune::BIW407::secondOrderContraction(cauchyStressInkrement,
                                                     virtDeltStrain) *
                quadPoint.weight() * integrationElement;
          }
        }
        residualVector[dim * row + i] -=
            Dune::BIW407::secondOrderContraction(cauchyStresses,
                                                 realDeltStrain) *
            quadPoint.weight() * integrationElement;
      }
    }
  }

  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      std::cout << elementMatrix[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template <class Basis, class Matrix>
void setOccupationPattern(const Basis &basis, Matrix &matrix) {
  enum { dim = Basis::GridView::dimension };
  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  MatrixIndexSet nb;
  nb.resize(basis.size(), basis.size());
  // A view on the FE basis on a single element
  auto localView = basis.localView();
  // Loop over all leaf elements
  for (const auto &element : elements(basis.gridView())) {
    // Bind the local view to the current element
    localView.bind(element);
    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i = 0; i < localView.size(); i++) {
      // Global index of the i-th local degree of freedom of the current element
      auto row = localView.index(i);
      for (size_t j = 0; j < localView.size(); j++) {
        // Global index of the j-th local degree of freedom of the current
        // element
        auto col = localView.index(j);
        nb.add(row[0], col[0]); // im ersten Multiindex steht der block index!
      }
    }
  }
  // Give the matrix the occupation pattern we want.
  using namespace Indices;
  nb.exportIdx(matrix);
}

template <class Matrix, class MultiIndex>
decltype(auto) matrixEntry(Matrix &matrix, const MultiIndex &row,
                           const MultiIndex &col) {
  using namespace Indices;
  return matrix[row[0]][col[0]][row[1]][col[1]];
}

template <class Vector, class MultiIndex>
decltype(auto) vectorEntry(Vector &vector, const MultiIndex &row) {
  using namespace Indices;
  return vector[row[0]][row[1]];
}

// Assemble the Laplace stiffness matrix on the given grid view
template <class Matrix, class Vector, class CurvedGridView,
          class DisplacementGridFunction>
void assembleStiffnessMatrix(const CurvedGridView &curvedGridView,
                             const DisplacementGridFunction &displacements,
                             Matrix &matrix, Vector &rhs) {
  // Set all entries to zero
  matrix = 0;
  rhs = 0;
  // A view of the basis
  auto localBasisView = displacements.basis().localView();
  // A view of the Displacement Function
  auto localDisplacementFunction = localFunction(displacements);
  auto material = Dune::BIW407::StVenantKirchhoffMaterial<2>(100.0, 170.0);

  // Traverse the curved grid
  for (const auto &element : elements(curvedGridView)) {
    // std::cout << element.type() << std::endl;
    //  Bind the local FE basis view to the current element, as well as the
    //  local Function
    localBasisView.bind(element);
    localDisplacementFunction.bind(element);

    // Now let’s get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    // displacement function must derivative from
    Dune::Matrix<double> elementMatrix;
    std::vector<double> elementResidualVector;
    assembleElementStiffnessMatrix(localBasisView, localDisplacementFunction,
                                   material, elementMatrix,
                                   elementResidualVector);
    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i = 0; i < elementMatrix.N(); i++) {
      // The global index of the i-th local degree of freedom
      // of the current element
      auto row = localBasisView.index(i);
      for (size_t j = 0; j < elementMatrix.M(); j++) {
        auto col = localBasisView.index(j);
        matrixEntry(matrix, row, col) += elementMatrix[i][j];
      }
      vectorEntry(rhs, row) += elementResidualVector[i];
    }
  }
}

int main(int argc, char **argv) {
  // Maybe initialize MPI
  MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);
  // World Dimension
  constexpr int dim = 2;
  // Linear Ansatz Order
  constexpr int p = 1;

  /////////////////////////////////////////////////////////
  // GRID DATA
  /////////////////////////////////////////////////////////
  // With this the Grid is Read and the geometry is setup.
  // Unfortunately I need to use a _moving Grid_
  // But the reference Grid is still a YaspGrid for startes
  // We use YaspGrid for starters
  using RefGrid = Dune::YaspGrid<2>;
  RefGrid refGrid{{1.0, 1.0}, {2, 2}};

  auto refGridView = refGrid.leafGridView();
  using RefGridView = decltype(refGridView);

  //////////////////////////////////////////////////////////
  // Function Space Basis
  //////////////////////////////////////////////////////////

  // this is not a GlobalDiscreteGridFuntion, only a GridViewable Function,
  auto spatialCoordinates =
      discreteGridViewFunction<2>(refGrid.leafGridView(), p);

  // Interpolate the initial Geometry
  auto unityFunctor = [](const auto &u) { return u; };
  Functions::interpolate(spatialCoordinates.basis(),
                         spatialCoordinates.coefficients(), unityFunctor);

  //////////////////////////////////////////////////
  // Moving Mesh setup
  //////////////////////////////////////////////////

  // Now we have the Actual Grid, on which we do our calculations
  // The MeshDisplacements are given via a std::ref to circumvent a copy.
  CurvedGrid grid{refGrid, std::ref(spatialCoordinates)};

  // We now contruct the gridView of the CurvedGrid
  auto gridView = grid.leafGridView();
  using GridView = decltype(gridView);

  using Coordinate = GridView::Codim<0>::Geometry::GlobalCoordinate;

  // constructing the basis:
  // which is a function on the reference geometry
  using namespace Functions::BasisFactory;
  auto lagrangeBasis =
      makeBasis(gridView, power<dim>(lagrange<p>(), blockedInterleaved()));
  // Blocked interleaved allows for efficent sparsity structure.

  ///////////////////////////////////////////////////////////
  // Setting up Storage and Matrizes
  ///////////////////////////////////////////////////////////

  // Every Node has `dim` Degrees of Freedom.
  using DisplacementRange = FieldVector<double, dim>;
  // Nodal Matrix
  using DisplacementMatrix = FieldMatrix<double, dim, dim>;
  // Multiindex: Nodes x Num Dimension
  using DisplacementVector = BlockVector<DisplacementRange>;
  // Double the Multiindex: One Block for each Node of Dimension (dim x dim)
  //
  using Matrix = BCRSMatrix<DisplacementMatrix>;
  // Marking displacement
  using DisplacementBitVector = std::vector<std::array<char, dim>>;

  // Declare the Boundary Vector, and the Right-Hand Side
  DisplacementBitVector isBoundary;
  DisplacementVector rhs;
  Matrix stiffnessMatrix;

  // 'rhsBackend' now connects rhs with `dune-functions`
  auto rhsBackend = Functions::istlVectorBackend(rhs);
  rhsBackend.resize(lagrangeBasis);

  rhs = 0;
  DisplacementVector u = rhs;
  DisplacementVector uIncrement = rhs; // Copy

  auto isBoundaryBackend = Functions::istlVectorBackend(isBoundary);
  isBoundaryBackend.resize(lagrangeBasis);

  //// THIS is the function that lives on the current geometry and will serve as
  ///the displacements from the initial geometry
  auto displacementFunction =
      Functions::makeDiscreteGlobalBasisFunction<DisplacementRange>(
          lagrangeBasis, u);

  // This only works if the basis from the geometry matches with the basis of
  // the Displacements
  auto displacementFunctionMaterial =
      Functions::makeDiscreteGlobalBasisFunction<DisplacementRange>(
          spatialCoordinates.basis(), u);
  //////////////

  using namespace Indices;
  for (auto &&b0i : isBoundary)
    for (std::size_t j = 0; j < b0i.size(); ++j)
      b0i[j] = false;

  // Declare the type of the stiffness matrix

  // Set matrix size and occupation pattern
  setOccupationPattern(displacementFunction.basis(), stiffnessMatrix);

  ////////////////////////////////////////////////
  // Even more Setup
  DisplacementVector initialX = spatialCoordinates.coefficients(); // Copy

  // std::cout << "Assembling step." << std::endl;
  // Funktion die den Dirichlet-Rand vorgibt
  auto &&g = [](Coordinate xmat) {
    return DisplacementRange{
        -0.0, -0.01}; //(std::abs(xmat[1] - 1.0) < 1e-8) ? 0.00002 : 0.0};
  };

  // Convenience Function for Boundary DOFs
  // Currently Marks all Boundary dofs
  Functions::forEachBoundaryDOF(displacementFunction.basis(),
                                [&](auto &&index) {
                                  // vertices(gridView)
                                  isBoundaryBackend[index] = true;
                                });

  // Start the NEWTON LOOP
  int iter_num = 0;
  do {

    uIncrement = 0;

    iter_num++;

    // modify rhs to incorporate dirichlet values
    Functions::interpolate(displacementFunction.basis(), u, g, isBoundary);

    assembleStiffnessMatrix(gridView, displacementFunction, stiffnessMatrix,
                            rhs);

    // std::cout << "Dirichlet step." << std::endl;

    // stiffnessMatrix.mmv(u,rhs);

    // std::cout << stiffnessMatrix. << std::endl;

    // Modify Stiffness Matrix to incorporate Dirichletvalues
    {
      auto localView = displacementFunction.basis().localView();
      for (const auto &element : elements(gridView)) {
        localView.bind(element);
        for (size_t i = 0; i < localView.size(); ++i) {
          auto row = localView.index(i);
          // If row corresponds to a boundary entry,
          // modify it to be an identity matrix row.
          if (isBoundaryBackend[row]) {
            for (size_t j = 0; j < localView.size(); ++j) {
              auto col = localView.index(j);
              matrixEntry(stiffnessMatrix, row, col) = (i == j) ? 1.0 : 0;
            }
            vectorEntry(rhs, row) = vectorEntry(u, row);
          }
        }
      }
    }

    std::cout << rhs << std::endl;

    //////////////////////////////////////////////////
    // Solving the ~linear~ system
    //////////////////////////////////////////////////

    // Turn the matrix into a linear operator
    // std::cout << "Solving step." << std::endl;

    MatrixAdapter<Matrix, DisplacementVector, DisplacementVector>
        stiffnessOperator(stiffnessMatrix);

    // Sequential incomplete LU decomposition as the preconditioner
    SeqILU<Matrix, DisplacementVector, DisplacementVector> preconditioner(
        stiffnessMatrix, 1.0); // Relaxation factor
                               // Preconditioned conjugate gradient solver
    CGSolver<DisplacementVector> cg(stiffnessOperator, preconditioner,
                                    1e-13, // Desired residual reduction factor
                                    100,   // Maximum number of iterations
                                    2);    // Verbosity of the solver
    // Object storing some statistics about the solving process
    InverseOperatorResult statistics;
    // Solve!
    cg.apply(uIncrement, rhs, statistics);

    // Changing the Inkrement to exclude the dirichlet Values
    auto localView = displacementFunction.basis().localView();
    for (const auto &element : elements(gridView)) {
      localView.bind(element);
      for (size_t i = 0; i < localView.size(); ++i) {
        auto row = localView.index(i);
        // If row corresponds to a boundary entry,
        if (isBoundaryBackend[row]) {
          vectorEntry(uIncrement, row) = 0;
          // vectorEntry(spatialCoordinates.coefficients(), row)
          // +=vectorEntry(uIncrement, row);
        }
      }
    }
    u += uIncrement;

    spatialCoordinates.coefficients() = 0;
    spatialCoordinates.coefficients() += initialX;
    spatialCoordinates.coefficients() += u;

    // Interpolate new geometry

    std::cout << rhs.two_norm() << std::endl;
    std::cout << uIncrement << std::endl;
    std::cout << u << std::endl;
    std::cout << spatialCoordinates.coefficients() << std::endl;

    // Because the Displacements live on the reference configuration, the
    // reference _grid_ will be written.
    using DataCollector = Vtk::LagrangeDataCollector<RefGridView, p>;
    using Writer = VtkUnstructuredGridWriter<RefGridView, DataCollector>;
    Writer vtkWriter(refGridView);
    vtkWriter.addPointData(displacementFunctionMaterial, "displacement");
    vtkWriter.write("displacement-result.vtu");

  } while (rhs.two_norm() > 1e-10 && iter_num < 4);

  // std::cout << xIncrement << std::endl;
}
