
#include "dune/common/exceptions.hh"
#include <dune/biw4-07/biw4-07.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>
#include <ostream>
#include <random>

typedef std::mt19937
    MyRNG; // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val = 123456789; // populate somehow

MyRNG rng; // e.g. keep one global instance (per thread)

typedef Dune::FieldMatrix<double, 2, 2> M2x2;

void initialize() { rng.seed(seed_val); }

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

Dune::FieldMatrix<double, 2, 2> createDeformationGradient(double shear,
                                                          double elongation) {

  Dune::FieldMatrix<double, 2, 2> deformationGradient(0);
  deformationGradient[0][0] = 1.0 * elongation;
  deformationGradient[0][1] = shear;
  deformationGradient[1][0] = 0.0;
  deformationGradient[1][1] = 1.0;

  return deformationGradient;
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

double difference(const Dune::FieldMatrix<double, 2, 2> &A,
                  const Dune::FieldMatrix<double, 2, 2> &B) {
  return std::abs(A[0][0] - B[0][0]) + std::abs(A[0][1] - B[0][1]) +
         std::abs(A[1][0] - B[1][0]) + std::abs(A[1][1] - B[1][1]);
}

int main(int argc, char **argv) {
  initialize();

  std::uniform_real_distribution<double> distribution(-1.0, 1.0);

  auto material = Dune::BIW407::NeoHookeMaterial<2>(30.0, 70.0);

  for (int i = 0; i < 10000; i++) {

    double shear = distribution(rng);
    double elongation = distribution(rng);

    // Hello there
    const auto deformation = createDeformationGradient(shear, elongation);
    auto jacobian = deformation.determinant();

    M2x2 cauchyStresses(0);

    material.cauchyStresses(deformation, cauchyStresses);

    if (jacobian > 0) {
      for (int j = 0; j < 100; j++) {
        double rotation = distribution(rng);
        auto rotMatrix = createRotationMatrix(rotation);

        M2x2 sigmaStresses(0);
        auto defoRot = rotMatrix * deformation;
        material.cauchyStresses(defoRot, sigmaStresses);

        sigmaStresses = rotMatrix.transposed() * sigmaStresses * rotMatrix;

        if (difference(sigmaStresses, cauchyStresses) > 1e-5) {
          std::cerr << "Material is not Observer Invariant!" << std::endl;

          std::cerr << "Expected Stresses" << std::endl
                    << cauchyStresses << std::endl
                    << cauchyStresses.determinant() << std::endl;
          std::cerr << "Actual Stresses" << std::endl
                    << sigmaStresses << std::endl
                    << sigmaStresses.determinant() << std::endl;

          std::cerr << "Deformation Gradient" << std::endl
                    << deformation << std::endl
                    << deformation.determinant() << std::endl;
          std::cerr << "Rotated Deformation Gradient" << std::endl
                    << defoRot << std::endl
                    << defoRot.determinant() << std::endl;
          std::cerr << "Rotation Matrix" << std::endl << rotMatrix << std::endl;

          DUNE_THROW(Dune::Exception, "Material is not Observer Invariant");
        }
        // std::cout << jacobian << std::endl;
      }
    }
    // Material to test
  }
}