//
// Created by mitja on 4.6.2019.
//

#ifndef VAJA_VI_1_HELIUM_H
#define VAJA_VI_1_HELIUM_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

inline void operator/=(std::vector<double>& v, double a)
{
    for (double& x : v) x /= a;
}

inline std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2)
{
    std::vector<double> w(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), w.begin(), std::minus<>());
    return w;
}

inline double inner(const std::vector<double>& v1, const std::vector<double>& v2, double h = 1.0)
{
    return h * std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

inline double norm(const std::vector<double>& v, double h)
{
    return std::sqrt(inner(v, v, h));
}

class Helium
{
public:
    Helium(const std::vector<double>& u, double R);
    ~Helium();
    inline double epsilon() const {return m_e;}
    inline double norm() const {return ::norm(m_u, m_h);}
    std::vector<double> electrostatic();
    double energy();
    std::vector<double> psi();

private:
    inline void constructSpace(){ for (size_t n = 0; n < m_N; ++n) m_r[n] = (n + 1) * m_h; }
    void constructVarBasis();
    void constructHamiltonian();
    void updateDensity();
    std::vector<double> correlatic();
    void LDA_DFT();

    gsl_matrix* m_H;
    gsl_matrix* m_eigenVec;
    gsl_vector* m_eigenVal;
    gsl_eigen_symmv_workspace* m_eigenWorkspace;
    gsl_vector* m_eigenColumn;

    std::vector<std::vector<double> > m_basis;
    std::vector<double> m_u;
    std::vector<double> m_r;
    std::vector<double> m_V;
    std::vector<double> m_f;

    double m_h;
    double m_e;
    size_t m_N;

};


#endif //VAJA_VI_1_HELIUM_H
