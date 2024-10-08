#pragma once
#ifndef BIW4_07_HH
#define BIW4_07_HH
#include <cmath>

#include <dune/common/fvector.hh> 
#include <dune/common/fmatrix.hh>
#include <dune/biw4-07/currentConfigMaterial.hh>

// add your classes here
namespace Dune {
    namespace BIW407 {
        
        /** \brief Compute linearised strain at the identity from a given displacement gradient.
         *
         *  \param grad The gradient of the direction in which the linearisation is computed.
         *  \param strain The tensor to store the strain in.
         */
        template <int dim, class field_type=double>
        static void linearisedStrain(const Dune::FieldMatrix<field_type, dim, dim>& grad, Dune::FieldMatrix<double,dim,dim>& strain) {
            for (int i=0; i<dim ; ++i)
            {
                //strain[i][i] = grad[i][i];
                for (int j=0; j<dim; ++j)
                    strain[i][j] = 0.5*(grad[i][j] + grad[j][i]);
            }
        }

        template<int dim, class field_type=double>
        static double secondOrderContraction(const Dune::FieldMatrix<field_type, dim, dim>& first, const Dune::FieldMatrix<field_type, dim, dim>& second){
            double result = 0.0;
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    result += first[i][j] * second[i][j];
                }
            }
            return result;
        }


        template<int dim, class field_type=double>
        static Dune::FieldMatrix<field_type, dim, dim> inverse(const Dune::FieldMatrix<field_type, dim, dim>& A) {
            Dune::FieldMatrix<field_type, dim, dim> result(0);

            auto det = A.determinant();

            if (dim == 2) {
                result[0][0] = 1.0/det * A[1][1];
                result[0][1] = -1.0/det * A[1][0];
                result[1][0] = -1.0/det * A[0][1];
                result[1][1] = 1.0/det * A[0][0];
                return result;
            }
        }

        //! C = F * F^T
        template<int dim, class field_type=double>
        static void leftCauchyGreenStretch(const Dune::FieldMatrix<field_type, dim, dim>& deformationGradient,Dune::FieldMatrix<field_type, dim, dim>& leftCauchyGreen){
            for (int i=0; i<dim ; ++i) {
                for (int j=0; j<dim; ++j){
                    for (int k = 0; k < dim; ++k) {
                        leftCauchyGreen[i][j] += deformationGradient[i][k] * deformationGradient[j][k];
                    }
                }
            }
        }

        //! b = F^T * F
        template<int dim, class field_type=double>
        static void rightCauchyGreenStretch(const Dune::FieldMatrix<field_type, dim, dim>& deformationGradient,Dune::FieldMatrix<field_type, dim, dim>& rightCauchyGreen){
            for (int i=0; i<dim ; ++i) {
                for (int j=0; j<dim; ++j){
                    for (int k = 0; k < dim; ++k) {
                        rightCauchyGreen[i][j] += deformationGradient[k][i] * deformationGradient[k][j];
                    }
                }
            }
        }

        template<int dim>
        class NeoHookeMaterial : CurrentConfigMaterial<dim> {

            double _shearModulus, _bulkModulus;

            /*
                Constructing the Material from the given values
            */
            public:
                NeoHookeMaterial(double shearModulus, double bulkModulus) :
                    _shearModulus(shearModulus),
                    _bulkModulus(bulkModulus) {};

            double strainEnergyDensity(const Dune::FieldMatrix<double, dim, dim>& deformationGradient) const {
                double jacobian = deformationGradient.determinant();

                Dune::FieldMatrix<double, dim, dim> rightCauchyGreen(0.0);

                rightCauchyGreenStretch(deformationGradient, rightCauchyGreen);
                double traceC = 0.0;
                for (int i = 0; i < dim; i++) {
                    traceC += rightCauchyGreen[i][i];
                }
                traceC += 1.0; // Plane Strain


                return _shearModulus * 0.5 * ( std::pow(jacobian, -2.0/3.0) * traceC - 3.0 ) 
                + _bulkModulus * 0.25 * (jacobian*jacobian - 1.0 - 2.0 * std::log(jacobian));
            }

            void secondPKStresses(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, Dune::FieldMatrix<double, dim, dim>& secondPKStress) const {
                Dune::FieldMatrix<double, dim, dim> inverseDeformation = deformationGradient;
                inverseDeformation.invert();

                auto jacobian = deformationGradient.determinant();

                auto inverseRightCauchyGreen = inverseDeformation * inverseDeformation.transposed();
                auto rightCauchyGreen = deformationGradient.transposed() * deformationGradient;

                double traceC = 0.0;
                for (int i = 0; i<dim; i++) {
                    traceC += rightCauchyGreen[i][i];
                }
                traceC += 1.0;

                secondPKStress = 0.0;

                for (int i = 0; i <dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        secondPKStress[i][j] += (
                        _bulkModulus*0.5*(jacobian*jacobian - 1.0)*inverseRightCauchyGreen[i][j]
                        -_shearModulus * std::pow(jacobian,-2.0/3.0) * 1.0/3.0 * inverseRightCauchyGreen[i][j]*traceC
                        );
                    }
                    secondPKStress[i][i] +=_shearModulus * std::pow(jacobian,-2.0/3.0);
                }
            }

            /**
                \brief Computes the Cauchy Stresses from a given deformation gradient.

                T

                @param deformationGradient The given Deformation Gradient
                @param cauchyStress The FieldMatrix, in which the cauchy stress should be stored in
            */
            void cauchyStresses(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, Dune::FieldMatrix<double, dim, dim>& cauchyStress) const {
                double jacobian = deformationGradient.determinant();
                
                Dune::FieldMatrix<double, dim, dim> leftCauchy(0);
                double traceb = 0.0;
                
                leftCauchyGreenStretch(deformationGradient, leftCauchy);
                for (int i = 0; i < dim; i++) {
                    traceb += leftCauchy[i][i];
                }
                traceb += 1.0; // PLANE STRAIN CORRECTION!

                //std::cout << jacobian << std::endl;
                //std::cout << leftCauchy.determinant() << std::endl; 
                //std::cout << traceb << std::endl; 


                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++){
                        if (i == j) {
                            cauchyStress[i][j] = 
                            this->_shearModulus*std::pow(jacobian,-5.0/3.0)*(leftCauchy[i][j] - 1.0/3.0 * traceb)
                            + this->_bulkModulus*(jacobian*jacobian -1.0)/(2.0 * jacobian);
                        } else {
                            cauchyStress[i][j] = this->_shearModulus*std::pow(jacobian,-5.0/3.0)*(leftCauchy[i][j]);
                        }
                    }
                }
            }

            /**
                \brief Returns the Cauchy Stress Inkrement, given the linearized incremental strains

                @param deformationGradient the current deformation gradient
                @param symGradient the linerized incremental strains in the current configuration
                @param cauchyStressInkrement the current increment
             */
            void cauchyStressInkrement(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, const Dune::FieldMatrix<double, dim, dim>& symGradient, Dune::FieldMatrix<double, dim, dim>& cauchyStressInkrement) const {
                double jacobian = deformationGradient.determinant();
                
                Dune::FieldMatrix<double, dim, dim> leftCauchy(0);
                double traceb = 0.0; // tr b
                double traceStrain = 0.0; // tr delt e = 1 : delt e

                leftCauchyGreenStretch(deformationGradient, leftCauchy);
                for (int i = 0; i < dim; i++) {
                    traceb += leftCauchy[i][i];
                    traceStrain += symGradient[i][i];
                }
                traceb += 1.0;// Plane Strain

                // All others depending on b
                double froubeniusProduct = secondOrderContraction(leftCauchy, symGradient); // b : delt e

                //std::cout << froubeniusProduct << std::endl;
                //std::cout << symGradient << std::endl;
                Dune::FieldMatrix<double, dim, dim> indicator = symGradient;

                for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++) {
                        indicator[i][j] = std::abs(symGradient[i][j]) > 1e-6 ? 1.0 : 0.0;
                    }
                

                // Populating the matrix
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        if (i== j) {
                            cauchyStressInkrement[i][j] =
                            - 2.0*this->_shearModulus/3.0 * std::pow(jacobian,-5.0/3.0) 
                            *(
                                ( leftCauchy[i][j] - 1.0/3.0 * traceb ) * traceStrain
                                - traceb * symGradient[i][j]
                                + froubeniusProduct
                            )
                            + this->_bulkModulus/(jacobian) * (
                                jacobian*jacobian * traceStrain
                                -(jacobian*jacobian - 1.0) * symGradient[i][j]
                            );
                        } else {
                            cauchyStressInkrement[i][j] =
                            -2.0* this->_shearModulus/3.0 * std::pow(jacobian,-5.0/3.0) 
                            *(
                                ( leftCauchy[i][j] ) * traceStrain
                                - traceb * symGradient[i][j]
                            )
                            + this->_bulkModulus/(jacobian) * (
                                -(jacobian*jacobian - 1.0) * symGradient[i][j]
                            );
                        }
                        cauchyStressInkrement[i][j] *= indicator[i][j];
                    }
                }
            }

        };

        template<int dim>
        class StVenantKirchhoffMaterial : CurrentConfigMaterial<dim>{

            double _shearModulus,_firstLameParameter;

            public:
                StVenantKirchhoffMaterial(double shearModulus, double firstLameParameter):
                    _shearModulus(shearModulus),
                    _firstLameParameter(firstLameParameter)
                {};

            double strainEnergyDensity(const Dune::FieldMatrix<double, dim, dim>& deformationGradient) const {
                
                Dune::FieldMatrix<double, dim, dim> rightCauchy(0);
                
                rightCauchyGreenStretch(deformationGradient, rightCauchy);
                // Euler Almansi
                for (int i = 0; i <dim; i++){
                    rightCauchy[i][i] -= 1.0;
                }

                rightCauchy *= 0.5;
                 

                double traceE = 0.0;
                double froebenius = 0.0;
                for (int i = 0; i <dim; i++){
                    traceE += rightCauchy[i][i];
                    for (int j = 0; j <dim; j++){
                        froebenius += rightCauchy[j][i] * rightCauchy[j][i];
                    }
                }
                return (
                    _firstLameParameter * 0.5 * traceE * traceE
                   + _shearModulus * froebenius
                );
            }

            void secondPKStresses(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, Dune::FieldMatrix<double, dim, dim>& PKIIStress) const {
                Dune::FieldMatrix<double, dim, dim> rightCauchy(0);
                
                PKIIStress = 0.0;

                rightCauchyGreenStretch(deformationGradient, rightCauchy);
                // Euler Almansi
                for (int i = 0; i <dim; i++){
                    rightCauchy[i][i] -= 1.0;
                }
                // Now its Euler Almansi Strains...
                rightCauchy *= 0.5;

                double traceE = 0.0;
                for (int i = 0; i <dim; i++){
                    traceE += rightCauchy[i][i];
                }

                for (int i = 0; i <dim; i++){
                    for (int j = 0; j <dim; j++){
                        PKIIStress[i][j] += 2.0*_shearModulus*rightCauchy[i][j];
                    }
                    PKIIStress[i][i] += _firstLameParameter * traceE;
                }


            }

            void cauchyStresses(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, Dune::FieldMatrix<double, dim, dim>& cauchyStress) const {
                double jacobian = deformationGradient.determinant();
                
                Dune::FieldMatrix<double, dim, dim> leftCauchy(0.0);
                Dune::FieldMatrix<double, dim, dim> leftCauchySquared(0.0);
                Dune::FieldMatrix<double, dim, dim> rightCauchyGreen(0.0);
                double traceb = 0.0;

                leftCauchyGreenStretch(deformationGradient, leftCauchy);
                
                leftCauchySquared = leftCauchy * leftCauchy;
                
                for (int i = 0; i < dim; i++) {
                    traceb += leftCauchy[i][i];
                }
                traceb += 1.0; // Plane Strain Correction!

                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        cauchyStress[i][j] = 1.0/jacobian *
                        (
                            (_firstLameParameter) * 0.5 * ( traceb - 3.0 ) * (leftCauchy[i][j])
                            + (_shearModulus) * ( (leftCauchySquared[i][j]) - (leftCauchy[i][j]))
                        );
                    }
                }
            }

            void cauchyStressInkrement(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, const Dune::FieldMatrix<double, dim, dim>& symGradient, Dune::FieldMatrix<double, dim, dim>& cauchyStressInkrement) const {
                double jacobian = deformationGradient.determinant();
                
                Dune::FieldMatrix<double, dim, dim> leftCauchy(0);
                Dune::FieldMatrix<double, dim, dim> transportCauchy(0);
                double froebenius = 0.0;

                
                leftCauchyGreenStretch(deformationGradient, leftCauchy);
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++){
                        froebenius += leftCauchy[i][j] * symGradient[i][j];
                    }
                }
                //traceE += 1.0; // plane strain

                // b * delt e * b
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        for (int k = 0; k < dim; k++) {
                            for (int l = 0; l < dim; l++) {
                                transportCauchy[i][j] += leftCauchy[i][k] * symGradient[k][l] * leftCauchy[l][j];
                            }
                        }
                    }
                }

                cauchyStressInkrement = 0;

                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        cauchyStressInkrement[i][j] = 1.0/jacobian *
                        (
                            _firstLameParameter *leftCauchy[i][j] * froebenius  
                            + 2.0*_shearModulus * ( transportCauchy[i][j])
                        );
                    }
                }


            }
            

        };


    }
}


#endif // BIW4_07_HH
