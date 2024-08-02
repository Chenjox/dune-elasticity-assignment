#ifndef BIW4_07_HH
#define BIW4_07_HH

#include <cmath>

#include <dune/common/fvector.hh> 
#include <dune/common/fmatrix.hh>


// add your classes here
namespace Dune {
    namespace BIW407 {
    
        /**
            \brief Computes the linear strain bmatrix from the shape function gradient
         */
        template <int dim, int symdim = (dim*(dim-1)/2), class field_type=double>
        void linearStrainBMatrix(const Dune::FieldVector<field_type, dim>& grad, const Dune::FieldMatrix<field_type, symdim, dim>& bMatrix) {
            bMatrix = 0;
            if (dim == 2) {
                bMatrix[0][0] = grad[0];
                bMatrix[1][1] = grad[1];
                bMatrix[2][0] = grad[1];
                bMatrix[2][1] = grad[0];
            }
        }

        template <int dim, int symdim = (dim*(dim-1)/2), class field_type=double>
        void linearStrainCMatrix(const Dune::FieldMatrix<field_type, symdim, symdim>& CMatrix, const field_type first_lame, const field_type second_lame) {
            CMatrix = 0;
            if (dim == 2){
                CMatrix[0][0] = 2.0*second_lame + first_lame;
                CMatrix[0][1] = first_lame;
                CMatrix[1][1] = 2.0*second_lame + first_lame;
                CMatrix[1][0] = first_lame;
                CMatrix[2][2] = 2.0*second_lame;
            }
        }

    }
}


#endif // BIW4_07_HH
