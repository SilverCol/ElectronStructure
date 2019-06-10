//
// Created by mitja on 4.6.2019.
//

#ifndef VAJA_VI_1_SPECTRALH_H
#define VAJA_VI_1_SPECTRALH_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "Helium.h"

class SpectralH: public Helium
{
public:
    SpectralH(const std::vector<double>& u, double R);
    ~SpectralH();

private:
    void constructVarBasis();
    void constructHamiltonian();
    void updateDensity() override;

    gsl_matrix* m_H;
    gsl_matrix* m_eigenVec;
    gsl_vector* m_eigenVal;
    gsl_eigen_symmv_workspace* m_eigenWorkspace;
    gsl_vector* m_eigenColumn;

    std::vector<std::vector<double> > m_basis;
};


#endif //VAJA_VI_1_SPECTRALH_H
