/*
 * utilities.h
 *
 *  Created on: Mar 7, 2020
 *      Author: cje
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <string>



/* Namespace containing:
 * Methods for approximating integrals
 *  - Extended Midpoint (Rectangle) Method
 *  - Trapezoidal Method
 *  - Simpsons Method
 *
 * All methods follow this structure:
 * @param lb: lower bound for the integral
 * @param ub: upper bound for the integral
 * @param N: number of iterations to run the algorithm
 * @param func: function to estimate integral for
 * @return estimated integral value
 */
namespace integral_est {

template <class T>
double rect(double lb, double ub, int N, T &func)
{
    // step_size:
    double h = (ub - lb) / N;

    // Initialise estimate
    double Ah = 0.0;

    // Accumulate estimate
    for (double mid = lb + 0.5 * h; mid < ub; mid+=h)
    {
        //std::cout << "mid x value: " << mid << std::endl;
        //std::cout << "A(h): " << h*func(mid) << std::endl;
        Ah += h*func(mid);
    }

    // Return with error estimate
    return Ah + 1.0 /(N*N);
}

template <class T>
double trpz(double lb, double ub,  int N, T &func)
{
    // step_size:
    double h = (ub - lb) / N;
    //std::cout << "step_size: " << h << std::endl;

    // Initialise the estimate
    double Ah = 0.5 * func(lb) * h;

    // Accumulate the estimate
    for (double mid = lb + h; mid < ub; mid+=h)
    {
        //std::cout << "mid: " << mid << std::endl;
        Ah += h * func(mid);
        //std::cout << "Ah: " << Ah << std::endl;
    }

    // Evaluate at upper bound
    //std::cout << "mid: " << ub << std::endl;
    Ah += 0.5 * func(ub) * h;
    //std::cout << "Ah: " << Ah << std::endl;

    // Return with error estimate
    return Ah + 1.0 /(N*N);
}

template <class T>
double simp(double lb, double ub,  int N, T &func)
{
    // step_size:
    double h = (ub - lb) / N;
    //std::cout << "step_size: " << h << std::endl;

    // Initialise the estimate
    double Ah = 1.0/3.0 * func(lb) * h;

    // Accumulate the estimate
    int i = 1;
    for (double mid = lb + h; mid < ub; mid+=h)
    {
        Ah += (2 + (i % 2) * 2)/3.0 * h * func(mid);
        i++;
    }

    // Evaluate at upper bound
    //std::cout << "mid: " << ub << std::endl;
    Ah += 1/3 * func(ub) * h;
    //std::cout << "Ah: " << Ah << std::endl;

    // Return with error estimate
    return Ah + 1.0/3.0 /(N*N*N*N);
}

}
/* End of namespace integral_est */




namespace util
{
/*
 * Load in dataset
 * @param data_name: name of the data file
 * @param x, y: data gets loaded into these vectors
 */
void loadDataset(const std::string &data_name, int param, VecDoub_O &x, VecDoub_O &y)
{
    ifstream data(data_name);
    for(int i = 0; i < param; i++) {
        data >> y[i];
        data >> x[i];
    }
}

void loadDataset(const std::string &data_name, int param, VecDoub_O &x, VecDoub_O &y, VecDoub_O &z, VecDoub_O &w)
{
    ifstream data(data_name);
    for(int i = 0; i < param; i++) {
        data >> x[i];
        data >> y[i];
        data >> z[i];
        data >> w[i];
    }
}


/*
 *  Calculate design matrix
 *  @param vector x: x-data for the current system
 *  @param basis_n: degree of the polynomial
 *  @return design Matrix A
 */
MatDoub getDesignMatrix(const VecDoub &x, int basis_n)
{
    // The design matrix
    MatDoub A(x.size(), basis_n);

    // Loop through the rows in the design matrix
    for (int i = 0; i < x.size(); i++)
    {
        // Loop through the columns in the design matrix
        for (int j = 0; j < basis_n; j++)
        {
            double basis = 1.0;
            // Use the basis functions to calculate elements of the design matrix
            for (int mul = 0; mul < j; mul++)
            {
                basis *= x[i];
            }
            A[i][j] = basis;
        }
    }
    return A;
}




/*
 * Gram-Schmidt method
 * @param input: list of linearly independent vectors that form a basis S
 * @return: list of orthonormal vectors which form an orthonormal basis for S
 */
std::vector<VecDoub> gramSchmidt(std::vector<VecDoub> &input)
{
    std::vector<VecDoub> orthonormals;
    orthonormals.push_back(input[0]/input[0].length());

    for (int i = 1; i < input.size(); i++)
    {
        VecDoub sum(input[0].size(), 0.0);

        for (int j = 0; j <= i - 1; j++)
        {
            double scal = orthonormals[j]*input[i];
            VecDoub vec = orthonormals[j]*scal;
            sum = sum + vec;
        }
        VecDoub e = input[i] - sum;
        if (e.length() > numeric_limits<Doub>::epsilon())
            orthonormals.push_back(e/e.length());
        else
            std::cout << "orthonormal vector had a length close to zero!" << std::endl;
    }
    return orthonormals;
}

/* Least Squares Solution x
 * Finds closest point in subpace to given point
 * @param orth_basis: orthonormal basis for the subspace
 * @param b: given point
 * @return b_LS: Least Square Solution
 */
VecDoub leastSquares(std::vector<VecDoub> &orth_basis, const VecDoub &b)
{
    VecDoub b_LS(5, 0.0);
    for (auto &u : orth_basis)
    {
        b_LS = b_LS + u * (u * b);
    }
    return b_LS;
}

/*
 * Calculate residual error
 * @param A: design matrix
 * @param b: right hand side b
 * @param x: x-coordinates of system
 * @return residual error
 */
double residualError(MatDoub &A, VecDoub &x, VecDoub &b)
{
    VecDoub tmp = A*x-b;
    return tmp.length()/b.length();
}

/*
 * Examine residuals with different thresholds
 * @param A: design matrix
 * @param b: right hand side b
 * @param coef: number of coefficients
 * @param thresholds: list of different thresholds to test for
 * @return none
 * @note prints out residual errors and solved coefficients for different values
 * of thresholds for comparison
 */
void examineResiduals(MatDoub &A, VecDoub &b, int coef, const std::vector<double> &thresholds)
{
    for (auto &trsh : thresholds)
    {
        std::cout << "-------------------------" << std::endl;
        std::cout << "SVD - trsh: " << trsh << std::endl;
        SVD SVDi(A);
        VecDoub aSVD(coef);
        SVDi.solve(b,aSVD, trsh);
        std::cout << "Solved Coefficients" << std::endl;
        aSVD.print();

        std::cout << "Singular values:" << std::endl;
        SVDi.w.print();

        std::cout << "Machine Precision: " << numeric_limits<Doub>::epsilon() <<std::endl;
        std::cout << "Reciprocal Condition number: " << SVDi.inv_condition() << std::endl;

        std::cout << "Residual error: " << residualError(A, aSVD, b) << std::endl << std::endl << std::endl;

    }
}

/*
 * Error estimates on the solution
 * @param A: design matrix A
 * @param thresh: threshold for the singular values
 * @return a list of errors estimates for each of the fitted coefficients
 */
VecDoub errorEstimates(const MatDoub &A, double thresh = -1.)
{
    SVD SVDerror(A);
    VecDoub w = SVDerror.w;
    MatDoub v = SVDerror.v;
    VecDoub res(v.ncols(), 0.0);
    double tsh = (thresh >= 0. ? thresh : SVDerror.tsh); // Set threshold to user specified or default for the SVD object

    for (int j = 0; j < v.ncols(); j++)
    {
        double sum = 0.0;
        for (int i = 0; i < w.size(); i++)
        {
            if (w[i] > tsh)
            {
                double div = v[j][i]/w[i];
                sum += div*div;
            }
        }
        res[j] = sqrt(sum);
    }

    return res;
}

/*
 * Table with Richardson Extrapolation
 */
template <typename T>//, typename P>//, class P>
void doSummaryTable(int max_iter, double order, T &func)//, P &integrator)
{
    int s = 16;
    // Make header
    std::cout << "i" << setw(s) << "A(hi)" << setw(s) << "A(hi-1)-A(hi)"
              << setw(s) << "Rich-alp^k" << setw(s) << "Rich-error"
              << setw(s) << "F-comp" << std::endl;

    // Initial step size
    double h1 = 1.0;
    double h2 = 1.0;
    double h3 = 1.0;

    // Results
    double Ah1 = 0.0;
    double Ah2 = 0.0;
    double Ah3 = 0.0;

    // bounds for the integrals
    double lb = 0.0;
    double ub = 1.0;

    int acc_f_comp = 0;

    for (int i = 1; i <= max_iter; ++i)
    {
        // Save previous results
        Ah3 = Ah2;
        Ah2 = Ah1;

        // Iterations to do to get Ah1:
        int N = (ub - lb) / h1;
        acc_f_comp += N;

        //std::cout << "N is: " << N << std::endl;
        Ah1 = integral_est::rect(lb, ub, N, func);

        std::cout << i << setw(s) << Ah1 << setw(s);

        if (i > 1)
            std::cout << (Ah2 - Ah1) << setw(s);
        else
            std::cout << "*" << setw(s);
        if (i > 2)
            //std::cout << ( (Ah3 - Ah2) / (Ah2 - Ah1) - pow(2.0, order) ) << setw(s);
            std::cout << ( (Ah3 - Ah2) / (Ah2 - Ah1) ) << setw(s);
        else
            std::cout << "*" << setw(s);

        std::cout << ( (Ah2 - Ah1) / (pow(2,order) - 1) ) << setw(s) << acc_f_comp << std::endl;


        h3 = h2;
        h2 = h1;
        h1 *= 0.5;
    }
}

}
/* End of namespace util */









namespace as1 {

MatDoub getDesignMatrix(const VecDoub &th1, const VecDoub &th2)
{
    int i, j;

    // The design matrix
    MatDoub A(th1.size()*2, 4, 0.0);

    // Loop through all the x rows in the design matrix
    for (i = 0; i < th1.size(); i++)
    {
        j = i*2;
        A[j][0] = 1.0;
        A[j][2] = cos(th1[i]);
        A[j][3] = cos(th1[i]+th2[i]);

        j++;
        A[j][1] = 1.0;
        A[j][2] = sin(th1[i]);
        A[j][3] = sin(th1[i]+th2[i]);
    }

    return A;
}

VecDoub getZ(const VecDoub &x, const VecDoub &y)
{
    int i;
    VecDoub z(x.size()*2, 0.0);
    for (i = 0; i < x.size(); i++)
    {
        z[i*2] = x[i];
        z[i*2+1] = y[i];
    }

    return z;
}

}
/* End of namespace as1 */


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

#endif /* UTILITIES_H_ */
