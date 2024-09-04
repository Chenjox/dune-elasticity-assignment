
#pragma once


#include <dune/common/fvector.hh> 
#include <dune/common/fmatrix.hh>

namespace Dune {
    namespace BIW407 {

        template<int dim>
        class CurrentConfigMaterial {

            virtual double strainEnergyDensity(const Dune::FieldMatrix<double, dim, dim>& deformationGradient) const;
        
            virtual void cauchyStresses(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, Dune::FieldMatrix<double, dim, dim>& cauchyStress) const;

            virtual void cauchyStressInkrement(const Dune::FieldMatrix<double, dim, dim>& deformationGradient, const Dune::FieldMatrix<double, dim, dim>& symGradient, Dune::FieldMatrix<double, dim, dim>& cauchyStressInkrement) const;
        };

    }
}