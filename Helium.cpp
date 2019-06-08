//
// Created by mitja on 4.6.2019.
//

#include <iostream>
#include <gsl/gsl_sf_laguerre.h>
#include "Helium.h"


namespace
{
    const size_t M = 100;
    const double epsilon = 1e-2;
    const double F = std::pow(3.141592653589793 / 6, 1.0/3.0);
    const double A = F * .0545 * std::pow(8.0 * 3.141592653589793 / 3.0, -1.0/3.0);
    const double B = 11.4 * std::pow(8.0 * 3.141592653589793 / 3.0, 1.0/3.0);
}

Helium::Helium(const std::vector<double>& u, double R) :
m_basis(M, std::vector<double>(u.size())),
m_u(u),
m_r(u.size()),
m_V(u.size()),
m_f(u.size()),
m_h(R/u.size()),
m_E(0.0),
m_N(u.size())
{
    m_H = gsl_matrix_calloc(M, M);
    m_eigenVec = gsl_matrix_alloc(M, M);
    m_eigenVal = gsl_vector_alloc(M);
    m_eigenWorkspace = gsl_eigen_symmv_alloc(M);
    m_eigenColumn = gsl_vector_alloc(M);
    constructSpace();
    constructVarBasis();
    updateDensity();
    //LDA_DFT();
}

Helium::~Helium()
{
    gsl_eigen_symmv_free(m_eigenWorkspace);
    gsl_matrix_free(m_H);
    gsl_matrix_free(m_eigenVec);
    gsl_vector_free(m_eigenVal);
    gsl_vector_free(m_eigenColumn);
}

void Helium::constructVarBasis()
{
    for (size_t n = 1; n < M + 1; ++n)
    {
        for (size_t r = 0; r < m_N; ++r)
        {
            m_basis[n][r] = std::sqrt(4/std::pow(n, 5)) * std::exp(-m_r[r] / n) * gsl_sf_laguerre_n(n, 1.0, 2*m_r[r]/n);
        }
    }
}

void Helium::constructHamiltonian()
{
    for (size_t j = 0; j < M; ++j)
    {
        for (size_t k = 0; k < M; ++k)
        {
            double* element = gsl_matrix_ptr(m_H, j, k);

            for (size_t r = 0; r < m_N; ++r)
            {
                *element += m_basis[j][r] * m_r[r] * m_basis[k][r];
            }
            *element *= m_h;

            if (j == k) *element -= 1/(2 * j*j);
        }
    }
}

void Helium::updateDensity()
{
    constructHamiltonian();
    
    gsl_eigen_symmv(m_H, m_eigenVal, m_eigenVec, m_eigenWorkspace);
    gsl_eigen_symmv_sort(m_eigenVal, m_eigenVec, GSL_EIGEN_SORT_VAL_ASC);

    gsl_matrix_get_col(m_eigenColumn, m_eigenVec, 0);
    std::transform(m_u.begin(), m_u.end(), m_u.begin(), [](double u){return 0.0;});
    for (size_t n = 0; n < M; ++n)
    {
        std::transform(m_u.begin(), m_u.end(), m_basis.begin(), m_u.begin(),
                [this, n](double u, std::vector<double>& b)
                {
                    return u + gsl_vector_get(m_eigenColumn, n) * b[n];
                });
    }

    m_E = gsl_vector_get(m_eigenVal, 0);
}

std::vector<double> Helium::electrostatic()
{
    if (std::abs(norm() - 1.0) > epsilon*epsilon)
    {
        throw std::runtime_error("State not normalised.");
    }

    for (size_t r = 0; r < m_u.size(); ++r)
    {
        m_f[r] = - m_h*m_u[r]*m_u[r] / (12*(r+1));
    }

    std::vector<double> U(m_u.size(), 0.0);
    U[1] = m_h + (7*m_f[0] + 6*m_f[1] - m_f[2]) / 2;

    for (size_t n = 2; n < U.size(); ++n)
    {
        U[n] = 2*U[n-1] - U[n-2] + m_f[n] + 10*m_f[n-1] + m_f[n-2];
    }

    double k = (1 - U.back()) / U.size();
    for (size_t r = 0; r < U.size(); ++r)
    {
        U[r] += k*(r + 1);
    }

    return U;
}

void Helium::LDA_DFT()
{
    while(true)
    {
        std::vector<double> comparison(m_u);
        std::vector<double> U = electrostatic();
        for (size_t n = 0; n < m_u.size(); ++n)
        {
            m_f[n] = std::pow(m_u[n] / ((n+1)*m_h), 2.0/3.0);
        }
        for (size_t n = 0; n < m_V.size(); ++n)
        {
            m_V[n] = 2 * (U[n] - 1) / ((n + 1) * m_h) - F * m_f[n] - A * std::log(1 + B * m_f[n]);
        }

        updateDensity();
        if (::norm(m_u - comparison, m_h) < epsilon) return;
    }
}
