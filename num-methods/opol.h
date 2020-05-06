#ifndef OPOL_H
#define OPOL_H

// Operator Overloads

/*
 * Operator/ for matrix and scalar
 * @param A: matrix A
 * @param s: scalar
 * @return matrix with division applied or error if scalar is too small
 */
MatDoub operator/(const MatDoub &A, const double &s)
{
    if (s < numeric_limits<Doub>::epsilon())
    {
        std::cerr << "in matrix divide: s is smaller the machine precision!" << std::endl;
    }
    double div = 1/s;

    MatDoub res( A.nrows(), A.ncols());
    for(int i=0;i < A.nrows(); i++)
    {
        for(int j=0; j < A.ncols(); j++)
        {
            res[i][j] = A[i][j]*div;
        }
    }
    return res;
}


/*
 * Operator+ for std::vector and VecDoub
 * @param stdv: std::vector
 * @param v: NR VecDoub
 * @return std::vector
 */
template <typename T>
std::vector<double> operator+(const std::vector<T> &stdv, const VecDoub &v)
{
    if(stdv.size() != v.size())
        std::cerr << "In std::vector + VecDoub: different size of vectors!" << std::endl;

    std::vector<double> res;
    for (int i = 0; i < stdv.size(); i++)
    {
        res.push_back(stdv[i]+v[i]);
    }
    return res;
}


#endif // OPOL_H
