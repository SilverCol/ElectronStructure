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
    const double C1 = .0545;
    const double C2 = 11.4 * std::pow(8 * 3.141592653589793 / 3, 1.0/3.0);
    const double C = std::pow(6 / 3.141592653589793, 1.0/3.0);
}

Helium::Helium(const std::vector<double>& u, double R) :
m_basis(M, std::vector<double>(u.size())),
m_u(u),
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
    //         std::cout << ::inner(b1, b2, m_h) << std::endl;
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
            m_basis[n - 1][r] = m_r[r] * std::exp(-m_r[r] / n)
                    * gsl_sf_laguerre_n(n - 1, 1.0, 2*m_r[r]/n);
        }
        m_basis[n - 1] /= ::norm(m_basis[n - 1], m_h);
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
                *element += m_basis[j][r] * m_V[r] * m_basis[k][r];
            }
            *element *= m_h;

            size_t n = j + 1;
            if (j == k) *element -= 1/(2 * std::pow(n, 2));
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
        std::transform(m_u.begin(), m_u.end(), m_basis[n].begin(), m_u.begin(),
                [this, n](double u, double b)
                {
                    return u + gsl_vector_get(m_eigenColumn, n) * b;
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

    for (size_t n = 0; n < m_N; ++n)
    {
        m_f[n] = - std::pow(m_h*m_u[n], 2) / (12 * m_r[n]);
    }

    std::vector<double> U(m_N, 0.0);
    U[1] = m_h + (7*m_f[0] + 6*m_f[1] - m_f[2]) / 2;

    for (size_t n = 2; n < m_N; ++n)
    {
        U[n] = 2*U[n-1] - U[n-2] + m_f[n] + 10*m_f[n-1] + m_f[n-2];
    }

    double k = (1 - U.back()) / m_r.back();
    for (size_t r = 0; r < m_N; ++r)
    {
        U[r] += k * m_r[r];
    }

    return U;
}

std::vector<double> Helium::correlatic()
{
    std::vector<double> v(m_N);
    std::transform(m_u.begin(), m_u.end(), m_r.begin(), v.begin(), [](double u, double r)
    {
        return - C * std::pow(std::abs(u / r), 2.0/3.0)
                + C1 * std::log(1 + C2 * std::pow(std::abs(u / r), 2.0/3.0));
    });

    return v;
}

void Helium::LDA_DFT()
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

double Helium::energy()
{
    std::vector<double> U = electrostatic();
    std::vector<double> V = correlatic();

    double energy = 0.0;
    for (size_t n = 0; n < m_N; ++n)
    {
        energy -= 2 * U[n] * std::pow(m_u[n], 2) / m_r[n];
        energy -= V[n] * std::pow(m_u[n], 2) / 2;
    }
    energy *= m_h;
    return 2*m_e + energy;
}

std::vector<double> Helium::psi()
{
    std::vector<double> psi(m_N);
    std::transform(m_u.begin(), m_u.end(), m_r.begin(), psi.begin(), std::multiplies<>());
    return psi;
}
