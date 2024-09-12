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
uint32_t seed_val = 876556; // populate somehow

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

Dune::FieldVector<double, 2>
calcInvariants(const Dune::FieldMatrix<double, 2, 2> &mat) {
  double det = mat.determinant();
  double trace = mat[0][0] + mat[1][1];

  Dune::FieldVector<double, 2> result = {det, trace};

  return result;
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

  auto material = Dune::BIW407::NeoHookeMaterial<2>(100.0, 170.0);

  int entries[4][2] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};

  double perturb = 1e-7;
  for (int _i = 0; _i < 1000; _i++) {

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

    const auto deformation = deformationStep1;
    const auto jacobian = deformationStep1.determinant();

    // Check if it is physical
    if (jacobian != 0 && jacobian > 1e-4) {
        // std::cout << jacobian << std::endl;
      Dune::FieldMatrix<double, 2, 2> rightCauchyGreen(0);
      // const auto inverseDeformation = invertDeformationGradient(deformation);
      Dune::BIW407::rightCauchyGreenStretch(deformation, rightCauchyGreen);
      const auto corrDeformation = rightCauchyToDeformation(rightCauchyGreen);
      const auto inverseDeformation =
          invertDeformationGradient(corrDeformation);

      Dune::FieldMatrix<double, 2, 2> SecondPK(0);
      Dune::FieldMatrix<double, 2, 2> cauchy(0);

      material.secondPKStresses(corrDeformation, SecondPK);
      material.cauchyStresses(corrDeformation, cauchy);

      auto sigma = 1.0/jacobian * corrDeformation * SecondPK * corrDeformation.transposed();

      for (int m = 0; m < 4; m++) {
        auto entryi = entries[m][0];
        auto entryj = entries[m][1];

        if (std::abs(cauchy[entryi][entryj] -sigma[entryi][entryj]) > 1e-4) {
            std::cout << "Actual:" << std::endl
                  << cauchy << std::endl;
            std::cout << "Actual Invariants" << std::endl
                  << calcInvariants(cauchy) << std::endl;

            std::cout << "Numerically estimated" << std::endl
                  << sigma << std::endl;
            
            std::cout << "Numerical Invariants" << std::endl
                  << calcInvariants(sigma) << std::endl;

            std::cout << "Deformation Gradient" << std::endl
                  << corrDeformation << std::endl;

            std::cout << "Jacobian" << std::endl
                  << jacobian << std::endl;

            
            DUNE_THROW(Dune::Exception, "Derivative does not match!");
        }
      }

    }
  }
}