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


    }
}


#endif // BIW4_07_HH
