#pragma once
#ifndef BIW4_07_HH
#define BIW4_07_HH
#include <cmath>

#include <dune/common/fvector.hh> 
#include <dune/common/fmatrix.hh>


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
                strain(i,i) = grad[i][i];
                for (int j=i+1; j<dim; ++j)
                    strain(i,j) = 0.5*(grad[i][j] + grad[j][i]);
            }
        }

        template<int dim>
        class CurrentConfigMaterial {

            double _shearModulus, _bulkModulus;

            /*
                Constructing the Material from the given values
            */
            public:
                CurrentConfigMaterial(double shearModulus, double bulkModulus) :
                    _shearModulus(shearModulus),
                    _bulkModulus(bulkModulus) {};

            /**
                \brief Computes the Cauchy Stresses from a given deformation gradient.

                T

                @param deformationGradient The given Deformation Gradient
                @param cauchyStress The FieldMatrix, in which the cauchy stress should be stored in
            */
            void cauchyStresses(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, Dune::FieldMatrix<double, dim, dim>& cauchyStress) {
                double jacobian = deformationGradient.determinant();
                
                Dune::FieldMatrix<double, dim, dim> leftCauchy(0);
                double traceb = 0.0;
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        for (int k = 0; j < dim; k++){
                            leftCauchy[i][j] = deformationGradient[i][k] + deformationGradient[j][k];
                        }
                    }
                    traceb += leftCauchy[i][i];
                }

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
            void cauchyStressInkrement(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, const Dune::FieldMatrix<double, dim, dim>& symGradient, Dune::FieldMatrix<double, dim, dim> cauchyStressInkrement) {
                double jacobian = deformationGradient.determinant();
                
                Dune::FieldMatrix<double, dim, dim> leftCauchy(0);
                double traceb = 0.0; // tr b
                double traceStrain = 0.0; // tr delt e = 1 : delt e
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        for (int k = 0; j < dim; k++){
                            leftCauchy[i][j] = deformationGradient[i][k] + deformationGradient[j][k];
                        }
                    }
                    traceb += leftCauchy[i][i];
                    traceStrain += symGradient[i][i];
                }
                // All others depending on b
                double froubeniusProduct = 0; // b : delt e
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        froubeniusProduct += leftCauchy[i][j] * symGradient[i][j];
                    }
                }

                // Populating the matrix
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        if (i== j) {
                            cauchyStressInkrement[i][j] =
                            - this->_shearModulus/3.0 * std::pow(jacobian,-5.0/3.0) 
                            *(
                                ( leftCauchy[i][j] - 1.0/3.0 * traceb ) * traceStrain
                                - traceb * symGradient[i][j]
                                + froubeniusProduct
                            )
                            + this->_bulkModulus/(2.0 * jacobian) * (
                                jacobian*jacobian * traceStrain
                                -(jacobian*jacobian - 1.0) * symGradient[i][j]
                            );
                        } else {
                            cauchyStressInkrement[i][j] =
                            - this->_shearModulus/3.0 * std::pow(jacobian,-5.0/3.0) 
                            *(
                                ( leftCauchy[i][j] ) * traceStrain
                                - traceb * symGradient[i][j]
                            )
                            + this->_bulkModulus/(2.0 * jacobian) * (
                                -(jacobian*jacobian - 1.0) * symGradient[i][j]
                            );
                        }
                    }
                }



            }

        };


    }
}


#endif // BIW4_07_HH
