#include "dune/common/exceptions.hh"
#include <dune/biw4-07/biw4-07.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>
#include <iostream>
#include <ostream>
#include <random>

typedef std::mt19937
    MyRNG; // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val = 123456789; // populate somehow

MyRNG rng; // e.g. keep one global instance (per thread)

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

Dune::FieldMatrix<double, 2, 2>
createDeformationGradient(double shear, double elongation, double rotation) {

  Dune::FieldMatrix<double, 2, 2> deformationGradient(0);

  double sinRotation = std::sin(rotation);
  double cosRotation = std::cos(rotation);

  deformationGradient[0][0] = 1.0 * elongation * cosRotation;
  deformationGradient[0][1] = shear * (-sinRotation);
  deformationGradient[1][0] = sinRotation;
  deformationGradient[1][1] = 1.0 * cosRotation;

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

  double shear = distribution(rng);
  double elongation = distribution(rng);
  double rotation = distribution(rng);

  const auto deformation =
      createDeformationGradient(shear, elongation, rotation);
  auto unitTensor = Dune::FieldMatrix<double, 2, 2>(0);
  unitTensor[0][0] = 1.0;
  unitTensor[1][1] = 1.0;
  auto jacobian = deformation.determinant();

  double perturb = 1e-12;
  if (jacobian > 0) {
    const auto inverseDeformation = invertDeformationGradient(deformation);

    // C = F^T F
    Dune::FieldMatrix<double, 2, 2> rightCauchyGreen(0);
    // b = F F^T
    Dune::FieldMatrix<double, 2, 2> leftCauchyGreen(0);

    // Populate
    Dune::BIW407::leftCauchyGreenStretch(deformation, leftCauchyGreen);
    Dune::BIW407::rightCauchyGreenStretch(deformation, rightCauchyGreen);

    // Pushforward...
    auto rightCauchyIdentity =
        inverseDeformation.transposed() * rightCauchyGreen * inverseDeformation;
    if (difference(rightCauchyIdentity, unitTensor) > 1e-5) {
      std::cerr << "Identity not Identity" << std::endl;
      std::cerr << "It's" << std::endl << rightCauchyIdentity << std::endl;
      std::cerr << "Expected" << std::endl << unitTensor << std::endl;

      DUNE_THROW(Dune::Exception, "Right Cauchy Green Implementation at fault");
    }

    auto leftCauchyIdentity =
        inverseDeformation * leftCauchyGreen * inverseDeformation.transposed();
    if (difference(leftCauchyIdentity, unitTensor) > 1e-5) {
      std::cerr << "Identity not Identity" << std::endl;
      std::cerr << "It's" << std::endl << leftCauchyIdentity << std::endl;
      std::cerr << "Expected" << std::endl << unitTensor << std::endl;

      DUNE_THROW(Dune::Exception, "Right Cauchy Green Implementation at fault");
    }
  }
}