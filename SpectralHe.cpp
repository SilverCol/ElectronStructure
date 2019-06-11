//
// Created by mitja on 4.6.2019.
//

#include <iostream>
#include <gsl/gsl_sf_laguerre.h>
#include "SpectralHe.h"

namespace
{
    const size_t M = 10;
    const double tolerance = 1e-4;
}

SpectralHe::SpectralHe(const std::vector<double>& u, double R) :
Helium(u, R),
m_basis(M, std::vector<double>(m_N))
{
    m_H = gsl_matrix_calloc(M, M);
    m_eigenVec = gsl_matrix_alloc(M, M);
    m_eigenVal = gsl_vector_alloc(M);
    m_eigenWorkspace = gsl_eigen_symmv_alloc(M);
    m_eigenColumn = gsl_vector_alloc(M);
    constructVarBasis();
    // for (auto& b1 : m_basis)
    // {
    //     for (auto& b2 : m_basis)
    //     {
    //         std::cout << ::inner(b1, b2, m_h) << std::endl;
    //     }
    // }
    // updateDensity();
    LDA_DFT();
}

SpectralHe::~SpectralHe()
{
    gsl_eigen_symmv_free(m_eigenWorkspace);
    gsl_matrix_free(m_H);
    gsl_matrix_free(m_eigenVec);
    gsl_vector_free(m_eigenVal);
    gsl_vector_free(m_eigenColumn);
}

void SpectralHe::constructVarBasis()
{
    std::cout << "Constructing variational basis." << std::endl;
    for (size_t n = 1; n < M + 1; ++n)
    {
        std::cout << n << '/' << M << '\r' << std::flush;
        for (size_t r = 0; r < m_N; ++r)
        {
            m_basis[n - 1][r] = m_r[r] * std::exp(-m_r[r] / n)
                    * gsl_sf_laguerre_n(n - 1, 1.0, 2*m_r[r]/n);
        }
        m_basis[n - 1] /= ::norm(m_basis[n - 1], m_h);
    }
    std::cout << std::endl;
}

void SpectralHe::constructHamiltonian()
{
    for (size_t j = 0; j < M; ++j)
    {
        for (size_t k = 0; k < M; ++k)
        {
            double* element = gsl_matrix_ptr(m_H, j, k);
            *element = 0.0;

            for (size_t r = 0; r < m_N; ++r)
            {
                *element += m_basis[j][r] * m_V[r] * m_basis[k][r];
            }
            *element *= m_h;

            size_t n = j + 1;
            if (j == k) *element -= 1/(2 * std::pow(n, 2));
        }
    }
}

void SpectralHe::updateDensity()
{
    constructHamiltonian();

    gsl_eigen_symmv(m_H, m_eigenVal, m_eigenVec, m_eigenWorkspace);
    gsl_eigen_symmv_sort(m_eigenVal, m_eigenVec, GSL_EIGEN_SORT_VAL_ASC);

    gsl_matrix_get_col(m_eigenColumn, m_eigenVec, 0);
    std::transform(m_u.begin(), m_u.end(), m_u.begin(), [](double u){return 0.0;});
    for (size_t n = 0; n < M; ++n)
    {
        std::transform(m_u.begin(), m_u.end(), m_basis[n].begin(), m_u.begin(),
                [this, n](double u, double b)
                {
                    return u + gsl_vector_get(m_eigenColumn, n) * b;
                });
    }

    m_e = gsl_vector_get(m_eigenVal, 0);
}

void SpectralHe::LDA_DFT()
{
    std::cout << "Running LDA DFT." << std::endl;
    while(true)
    {
        std::cout << "Now at " << m_e << '\r' << std::flush;

        std::vector<double> comparison(m_u);
        std::vector<double> U = electrostatic();
        std::vector<double> V = correlatic();

        for (size_t n = 0; n < m_N; ++n)
        {
            m_V[n] = (2*U[n] - 1) / m_r[n] + V[n];
        }

        updateDensity();
        if (::norm(m_u - comparison, m_h) < tolerance) return;
    }
}
