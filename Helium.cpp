//
// Created by mitja on 4.6.2019.
//

#include <iostream>
#include <gsl/gsl_sf_laguerre.h>
#include "Helium.h"

namespace
{
    const size_t M = 10;
    const double tolerance = 1e-10;
    const double F = std::pow(3 / 3.141592653589793, 1.0/3.0);
    const double A = .0545 * std::pow(3 / (2 * 3.141592653589793), 2.0/3.0);
    const double B = 11.4 * std::pow(4 * 3.141592653589793 / 3, 1.0/3.0);
    const double C = std::pow(3 / (2 * std::pow(3.141592653589793, 2)), 1.0/3.0);
}

Helium::Helium(const std::vector<double>& u, double R) :
m_basis(M, std::vector<double>(u.size())),
m_psi(u),
m_r(u.size()),
m_V(u.size()),
m_f(u.size()),
m_h(R/u.size()),
m_e(0.0),
m_N(u.size())
{
    m_H = gsl_matrix_calloc(M, M);
    m_eigenVec = gsl_matrix_alloc(M, M);
    m_eigenVal = gsl_vector_alloc(M);
    m_eigenWorkspace = gsl_eigen_symmv_alloc(M);
    m_eigenColumn = gsl_vector_alloc(M);
    constructSpace();
    constructVarBasis();
    // for (auto& b1 : m_basis)
    // {
    //     for (auto& b2 : m_basis)
    //     {
    //         std::cout << ::inner(b1, b2, m_r, m_h) << std::endl;
    //     }
    // }
    // updateDensity();
    LDA_DFT();
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
    std::cout << "Constructing variational basis." << std::endl;
    for (size_t n = 1; n < M + 1; ++n)
    {
        std::cout << n << '/' << M << '\r' << std::flush;
        for (size_t r = 0; r < m_N; ++r)
        {
            m_basis[n - 1][r] = std::sqrt(4/std::pow(n, 5)) * std::exp(-m_r[r] / n)
                    * gsl_sf_laguerre_n(n - 1, 1.0, 2*m_r[r]/n);
        }
        m_basis[n - 1] /= ::norm(m_basis[n - 1], m_r, m_h);
    }
    std::cout << std::endl;
}

void Helium::constructHamiltonian()
{
    for (size_t j = 0; j < M; ++j)
    {
        for (size_t k = 0; k < M; ++k)
        {
            double* element = gsl_matrix_ptr(m_H, j, k);
            *element = 0.0;

            for (size_t r = 0; r < m_N; ++r)
            {
                *element += m_basis[j][r] * m_V[r] * m_basis[k][r] * std::pow(m_r[r], 2);
            }
            *element *= m_h;

            if (j == k) *element -= 1/(2 * std::pow(j + 1, 2));
        }
    }
}

void Helium::updateDensity()
{
    constructHamiltonian();

    gsl_eigen_symmv(m_H, m_eigenVal, m_eigenVec, m_eigenWorkspace);
    gsl_eigen_symmv_sort(m_eigenVal, m_eigenVec, GSL_EIGEN_SORT_VAL_ASC);

    gsl_matrix_get_col(m_eigenColumn, m_eigenVec, 0);
    std::transform(m_psi.begin(), m_psi.end(), m_psi.begin(), [](double u){return 0.0;});
    for (size_t n = 0; n < M; ++n)
    {
        std::transform(m_psi.begin(), m_psi.end(), m_basis[n].begin(), m_psi.begin(),
                [this, n](double psi, double b)
                {
                    return psi + gsl_vector_get(m_eigenColumn, n) * b;
                });
    }

    m_e = gsl_vector_get(m_eigenVal, 0);
}

std::vector<double> Helium::electrostatic()
{
    if (std::abs(norm() - 1.0) > tolerance)
    {
        throw std::runtime_error("State not normalized.");
    }

    for (size_t r = 0; r < m_N; ++r)
    {
        m_f[r] = - m_r[r] * std::pow(m_h*m_psi[r], 2) / 12;
    }

    std::vector<double> U(m_psi.size(), 0.0);
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

std::vector<double> Helium::correlatic()
{
    std::transform(m_psi.begin(), m_psi.end(), m_f.begin(), [](double psi)
    {
        return std::pow(std::abs(psi), 2.0/3.0);
    });

    std::vector<double> v(m_N);
    std::transform(m_f.begin(), m_f.end(), v.begin(), [](double psi)
    {
        // return -F * psi - A * std::log(1 + B * psi);
        return C * psi;
    });

    return v;
}

void Helium::LDA_DFT()
{
    std::cout << "Running LDA DFT." << std::endl;
    while(true)
    {
        std::cout << "Now at " << m_e << '\r' << std::flush;

        std::vector<double> comparison(m_psi);
        std::vector<double> U = electrostatic();
        std::vector<double> V = correlatic();

        for (size_t n = 0; n < m_N; ++n)
        {
            m_V[n] = (2*U[n] - 1) / m_r[n] - V[n];
        }

        updateDensity();
        if (::norm(m_psi - comparison, m_r, m_h) < tolerance) return;
    }
}

double Helium::energy()
{
    std::vector<double> U = electrostatic();
    std::vector<double> V = correlatic();

    double energy = 0.0;
    for (size_t n = 0; n < m_N; ++n)
    {
        energy -= 2 * U[n] * std::pow(m_psi[n], 2) * m_r[n];
        energy += V[n] * std::pow(m_psi[n] * m_r[n], 2) / 2;
    }
    energy *= m_h;
    return 2*m_e + energy;
}
