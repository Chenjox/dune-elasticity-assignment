
#include "dune/common/exceptions.hh"
#include <dune/biw4-07/biw4-07.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <random>

typedef std::mt19937
    MyRNG; // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val = 98765; // populate somehow

MyRNG rng; // e.g. keep one global instance (per thread)

void initialize() { rng.seed(seed_val); }

Dune::FieldMatrix<double, 2, 2> createDeformationGradient(double shear,
                                                          double elongation) {

  Dune::FieldMatrix<double, 2, 2> deformationGradient(0);
  deformationGradient[0][0] = 1.0 * elongation;
  deformationGradient[0][1] = shear;
  deformationGradient[1][0] = 0.0;
  deformationGradient[1][1] = 1.0;

  return deformationGradient;
}

Dune::FieldMatrix<double, 2, 2> createRotationMatrix(double rotation) {

  Dune::FieldMatrix<double, 2, 2> rotationMatrix(0);

  double cosRot = std::cos(rotation);
  double sinRot = std::sin(rotation);

  rotationMatrix[0][0] = cosRot;
  rotationMatrix[0][1] = -sinRot;
  rotationMatrix[1][0] = sinRot;
  rotationMatrix[1][1] = cosRot;

  return rotationMatrix;
}

Dune::FieldMatrix<double, 2, 2> invertDeformationGradient(
    const Dune::FieldMatrix<double, 2, 2> &deformationGradient) {

  double d = deformationGradient.determinant();
  Dune::FieldMatrix<double, 2, 2> inverse(0);

  inverse[0][0] = 1.0 / d * deformationGradient[1][1];
  inverse[0][1] = -1.0 / d * deformationGradient[1][0];
  inverse[1][0] = -1.0 / d * deformationGradient[0][1];
  inverse[1][1] = 1.0 / d * deformationGradient[0][0];

  return inverse;
}

Dune::FieldMatrix<double, 2, 2>
rightCauchyToDeformation(const Dune::FieldMatrix<double, 2, 2> &rightCauchy) {

  Dune::FieldMatrix<double, 2, 2> deformationGradient(0);

  double deter = rightCauchy.determinant();
  double trace = rightCauchy[0][0] + rightCauchy[1][1];

  double s = std::sqrt(deter);
  double t = std::sqrt(trace + 2.0 * s);

  deformationGradient[0][0] = 1.0 / t * (rightCauchy[0][0] + s);
  deformationGradient[0][1] = 1.0 / t * rightCauchy[0][1];
  deformationGradient[1][0] = 1.0 / t * rightCauchy[1][0];
  deformationGradient[1][1] = 1.0 / t * (rightCauchy[1][1] + s);

  return deformationGradient;
}

int main(int argc, char **argv) {
  initialize();

  int constexpr dim = 2;

  std::uniform_real_distribution<double> distribution(-2.0, 2.0);

  auto material = Dune::BIW407::StVenantKirchhoffMaterial<2>(100.0, 170.0);

  int entries[4][2] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};

  double perturb = 1e-7;
  for (int _i = 0; _i < 2; _i++) {

    ////
    // Step 1: Create a viable Deformation Gradient
    double shear = distribution(rng);
    double elongation = distribution(rng);
    double rotation = distribution(rng);

    // Hello there
    const auto rotMatrix = createRotationMatrix(rotation);
    auto deformationStep1 =
        rotMatrix * createDeformationGradient(shear, elongation);
    auto jacobianStep1 = deformationStep1.determinant();

    // Pure elongation in 1 direction
    if (_i == 0) {
      deformationStep1[0][0] = 2.0;
      deformationStep1[0][1] = 0.0;
      deformationStep1[1][0] = 0.0;
      deformationStep1[1][1] = 0.5;
    }
    if (_i == 1) {
      deformationStep1[0][0] = elongation;
      deformationStep1[0][1] = 0.0;
      deformationStep1[1][0] = 0.0;
      deformationStep1[1][1] = 1.0 / elongation;
    }
    // Pure elongation in 2 direction
    if (_i == 2) {
      deformationStep1[0][0] = 0.0;
      deformationStep1[0][1] = elongation;
      deformationStep1[1][0] = -1.0 / elongation;
      deformationStep1[1][1] = 0.0;
    }
    // Pure Shear
    if (_i == 3) {
      deformationStep1[0][0] = 1.0 + shear;
      deformationStep1[0][1] = shear;
      deformationStep1[1][0] = 0.0;
      deformationStep1[1][1] = 1.0;
    }

    if (jacobianStep1 < 0) {
      deformationStep1 *= jacobianStep1;
    }

    const auto deformation = deformationStep1;
    const auto jacobian = deformationStep1.determinant();

    // Check if it is physical
    if (jacobian != 0 && std::abs(jacobian) > 1e-4) {
      // std::cout << jacobian << std::endl;
      Dune::FieldMatrix<double, 2, 2> rightCauchyGreen(0);
      // const auto inverseDeformation = invertDeformationGradient(deformation);
      Dune::BIW407::rightCauchyGreenStretch(deformation, rightCauchyGreen);
      const auto corrDeformation = rightCauchyToDeformation(rightCauchyGreen);

      // Now we have a current deformation gradient

      const auto inverseDeformation =
          invertDeformationGradient(corrDeformation);

      Dune::FieldMatrix<double, 2, 2> middleIIPK(0);

      material.secondPKStresses(corrDeformation, middleIIPK);

      // S = 2 psi/C \implies S_ij = 2 * p psi/ p C_ij
      // it is stored with non-sensible values
      double materialTensor[2][2][2][2];
      for (int m = 0; m < 4; m++) {
        int entryi = entries[m][0];
        int entryj = entries[m][1];

        Dune::FieldMatrix<double, 2, 2> cauchyGreenIplu = rightCauchyGreen;
        Dune::FieldMatrix<double, 2, 2> cauchyGreenImin = rightCauchyGreen;
        cauchyGreenImin[entryi][entryj] -= perturb;
        cauchyGreenIplu[entryi][entryj] += perturb;

        Dune::FieldMatrix<double, 2, 2> plusIIPK(0);
        Dune::FieldMatrix<double, 2, 2> minuIIPK(0);

        const auto plusDeformation = rightCauchyToDeformation(cauchyGreenIplu);
        const auto minuDeformation = rightCauchyToDeformation(cauchyGreenImin);

        material.secondPKStresses(plusDeformation, plusIIPK);
        material.secondPKStresses(minuDeformation, minuIIPK);

        auto forwardDifference = (plusIIPK - middleIIPK) / perturb;
        auto centralDifference = (plusIIPK - minuIIPK) / (2.0 * perturb);
        auto backwardDifference = (middleIIPK - minuIIPK) / perturb;

        for (int n = 0; n < 4; n++) {
          int entryI = entries[n][0];
          int entryJ = entries[n][1];

          materialTensor[entryi][entryj][entryI][entryJ] =
              2.0*(forwardDifference[entryI][entryJ] +
               centralDifference[entryI][entryJ] +
               backwardDifference[entryI][entryJ]) /
              3.0;
        }
      }
      // now the material Tensor is initialised.
      // We now compute the push-forward,
      double spatialTensor[dim][dim][dim][dim];
      // init
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          for (int k = 0; k < dim; k++) {
            for (int l = 0; l < dim; l++) {
              spatialTensor[i][j][k][l] = 0.0;
            }
          }
        }
      }
      // push forward
      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
          for (int k = 0; k < dim; k++)
            for (int l = 0; l < dim; l++)
              for (int a = 0; a < dim; a++)
                for (int b = 0; b < dim; b++)
                  for (int c = 0; c < dim; c++)
                    for (int d = 0; d < dim; d++) {
                      spatialTensor[a][b][c][d] +=
                          1.0 / jacobian *
                          corrDeformation[a][i] * // The Equation
                          corrDeformation[b][j] * // must not
                          corrDeformation[c][k] * // be split!
                          corrDeformation[d][l] * materialTensor[i][j][k][l];
                    }
      // We now have the spatial material tensor.
      // Now we can test the current function implementation
      Dune::FieldMatrix<double, dim, dim> stressInkrement(0);
      Dune::FieldMatrix<double, dim, dim> sigmaInkrement(0);

      Dune::FieldMatrix<double, dim, dim> symGradient(0);

      for (int _j = 0; _j < 100; _j++) {

        double first = distribution(rng);
        double second = distribution(rng);
        for (int i = 0; i < dim; i++)
          for (int j = 0; j < dim; j++) {
            double take = 0.0;
            if (i == j) {
              take = first;
            } else {
              take = second + first;
            }
            symGradient[i][j] += 0.5 * take;
            symGradient[j][i] += 0.5 * take;
          }

        material.cauchyStressInkrement(corrDeformation, symGradient,
                                       stressInkrement);

        for (int i = 0; i < dim; i++)
          for (int j = 0; j < dim; j++)
            for (int k = 0; k < dim; k++)
              for (int l = 0; l < dim; l++) {
                sigmaInkrement[i][j] +=
                    spatialTensor[i][j][k][l] * symGradient[k][l];
              }

        for (int i = 0; i < dim; i++)
          for (int j = 0; j < dim; j++) {
            if (std::abs(1.0 - sigmaInkrement[i][j] / stressInkrement[i][j]) >
                0.2) {

              Dune::FieldMatrix<double, dim, dim> relErrorMatrix(0);
              for (int k = 0; k < dim; k++)
                for (int l = 0; l < dim; l++) {
                  relErrorMatrix[k][l] =
                      sigmaInkrement[k][l] / stressInkrement[k][l];
                }

              std::cout << "Actual:" << std::endl
                        << stressInkrement << std::endl;
              std::cout << "Numerically estimated" << std::endl
                        << sigmaInkrement << std::endl;
              std::cout << "Multiple Matrix" << std::endl
                        << relErrorMatrix << std::endl;

              std::cout << "At Deformation Gradient" << std::endl
                        << corrDeformation << std::endl;
              std::cout << "With symmetric strains" << std::endl
                        << symGradient << std::endl;
              DUNE_THROW(Dune::Exception,
                         "Hessian does not match strain energy density!");
            }
          }
      }
    }
  }
  // Material to test
}
