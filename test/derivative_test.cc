
#include <dune/biw4-07/biw4-07.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>
#include <iomanip>
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

  std::uniform_real_distribution<double> distribution(0.2, 9.0);

  auto material = Dune::BIW407::NeoHookeMaterial<2>(100.0, 200.0);

  int entries[4][2] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};

  for (int i = 0; i < 10000; i++) {

    double shear = distribution(rng);
    double elongation = distribution(rng);
    double rotation = distribution(rng);

    // Hello there
    const auto rotMatrix = createRotationMatrix(rotation);
    const auto deformation =
        rotMatrix * createDeformationGradient(shear, elongation);
    auto jacobian = deformation.determinant();

    double perturb = 1e-10;
    if (jacobian > 0) {
      // std::cout << jacobian << std::endl;
      Dune::FieldMatrix<double, 2, 2> rightCauchyGreen(0);
      // const auto inverseDeformation = invertDeformationGradient(deformation);
      Dune::BIW407::rightCauchyGreenStretch(deformation, rightCauchyGreen);
      const auto corrDeformation = rightCauchyToDeformation(rightCauchyGreen);
      const auto inverseDeformation =
          invertDeformationGradient(corrDeformation);

      double middleEnergy = material.strainEnergyDensity(corrDeformation);
      Dune::FieldMatrix<double, 2, 2> cauchyStresses(0);
      Dune::FieldMatrix<double, 2, 2> SecondPK(0);

      material.cauchyStresses(corrDeformation, cauchyStresses);

      // sigma * F^-T
      SecondPK = jacobian * inverseDeformation *
                 (cauchyStresses * inverseDeformation.transposed());
      // F^-1 * sigma * F^-T
      // SecondPK.leftmultiply(inverseDeformation);

      // S = 2 psi/C \implies S_ij = 2 * p psi/ p C_ij
      for (int m = 0; m < 4; m++) {
        int entryi = entries[m][0];
        int entryj = entries[m][1];
        Dune::FieldMatrix<double, 2, 2> forward = rightCauchyGreen;
        Dune::FieldMatrix<double, 2, 2> backward = rightCauchyGreen;
        forward[entryi][entryj] += perturb;
        backward[entryi][entryj] -= perturb;
        const auto forwardDefo = rightCauchyToDeformation(forward);
        const auto backwardDefo = rightCauchyToDeformation(backward);
        double forwardEnergy = material.strainEnergyDensity(forwardDefo);
        double backwardEnergy = material.strainEnergyDensity(backwardDefo);

        double forwardDifference =
            2.0 * (forwardEnergy - middleEnergy) / (perturb);
        double backwardDifference =
            2.0 * (middleEnergy - backwardEnergy) / (perturb);
        double centralDifference = (forwardEnergy - backwardEnergy) / perturb;

        if (std::abs(SecondPK[entryi][entryj] - centralDifference) > 10.0 ||
            std::abs(SecondPK[entryi][entryj] - forwardDifference) > 10.0 ||
            std::abs(SecondPK[entryi][entryj] - backwardDifference) > 10.0) {
          // Error itself
          std::cerr << std::setprecision(9);
          std::cerr << "Error: Derivative does not match at [" << entryi << ", "
                    << entryj << "]!" << std::endl;
          std::cerr << "Forward/Central/Backward Difference: "
                    << forwardDifference << "/" << centralDifference << "/"
                    << backwardDifference << std::endl;
          std::cerr << "Derivative: " << SecondPK[entryi][entryj] << std::endl;
          std::cerr << "Multiples: "
                    << SecondPK[entryi][entryj] / forwardDifference
                    << std::endl;

          // Stresses in Full
          std::cerr << "2nd PK: " << std::endl << SecondPK << std::endl;
          std::cerr << "Cauchy: " << std::endl << cauchyStresses << std::endl;
          std::cerr << "F: " << std::endl << corrDeformation << std::endl;
          std::cerr << "F+eps: " << std::endl << forwardDefo << std::endl;
          std::cerr << "F-eps: " << std::endl << backwardDefo << std::endl;
          std::cerr << "F^-1" << std::endl << inverseDeformation << std::endl;
          std::cerr << "det F" << std::endl << jacobian << std::endl;

          // Right Cauchy Greens
          std::cerr << "C: " << std::endl << rightCauchyGreen << std::endl;
          std::cerr << "C+eps: " << std::endl << forward << std::endl;
          std::cerr << "C-eps: " << std::endl << backward << std::endl;

          // Energys
          std::cerr << "Actual Energy: " << std::endl
                    << middleEnergy << std::endl;
          std::cerr << "Actual Energy+eps: " << std::endl
                    << forwardEnergy << std::endl;
          std::cerr << "Actual Energy-eps: " << std::endl
                    << backwardEnergy << std::endl;

          DUNE_THROW(Dune::Exception, "Derivative does not match!");
        }
      }
    }
    // Material to test
  }
}