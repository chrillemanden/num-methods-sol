/*
 * lessons.h
 * Solutions to Numerical Methods lessons spring 2020
 * 6th semester in Bsc in Engineering (Robot Systems)
 * University of Southern Denmark
 *
 *  Created on: Mar 7, 2020
 *      Author: cje
 */


#ifndef LESSONS_H
#define LESSONS_H

#include "nr3.h"
#include "utilities.h"
#include "roots.h"
#include "qrdcmp.h"
#include "ludcmp.h"
#include "roots_multidim.h"


#define pi 3.14159265358979323846

/*
 * Try to solve Filip and Pont datasets with
 * LU and Cholesky decomposition
 */
void lesson2()
{
    VecDoub xFilip(82); VecDoub yFilip(82);
    util::loadDataset("FilipData.dat",82,xFilip,yFilip);

    VecDoub xPont(40); VecDoub yPont(40);
    util::loadDataset("PontData.dat",40,xPont,yPont);

    // your code

    MatDoub APont = util::getDesignMatrix(xPont, 3);
    MatDoub ATPont = APont.transpose();
    MatDoub CPont = ATPont * APont;
    VecDoub cPont = ATPont * yPont;


    // Make the decompositions
    std::cout << "--- Pontius dataset ---" << std::endl;
    std::cout << "Cholesky" << std::endl;
    Cholesky CHPont(CPont);
    VecDoub aCHPont(cPont.size());
    CHPont.solve(cPont, aCHPont);
    aCHPont.print();

    std::cout << "LU" << std::endl;
    LUdcmp LUPont(CPont);
    VecDoub aLUPont(cPont.size());
    LUPont.solve(cPont, aLUPont);
    aLUPont.print();


    // Filip
    MatDoub AFilip = util::getDesignMatrix(xFilip, 11);
    MatDoub ATFilip = AFilip.transpose();
    MatDoub CFilip = ATFilip * AFilip;
    VecDoub cFilip = ATFilip * yFilip;

    std::cout << std::endl << std::endl;
    // Make the decompositions
    std::cout << "--- Filip dataset ---" << std::endl;
    std::cout << "LU" << std::endl;
    LUdcmp LUFilip(CFilip);
    VecDoub aLUFilip(cFilip.size());
    LUFilip.solve(cFilip, aLUFilip);
    aLUFilip.print();

    /*
     * Cholesky will fail on the Filip dataset!
     * Therefore commented out
     */
//    std::cout << "Cholesky" << std::endl;
//    Cholesky CHFilip(CFilip);
//    VecDoub aCHFilip(cFilip.size());
//    CHFilip.solve(cFilip, aCHFilip);
//    util::print(aCHFilip);

}

/*
 * Performing the Singular Value Decomposition on the dataset
 * Filip where Cholensky decomposition failed
 * Using the Gram-Schmidt method to calculate an orthonormal basis
 * Find the closest point to the subspace that was previously found an orthonormal basis for
 */
void lesson3()
{

    VecDoub xFilip(82); VecDoub yFilip(82);
    util::loadDataset("FilipData.dat",82,xFilip,yFilip);
    MatDoub AFilip = util::getDesignMatrix(xFilip, 11);

    std::cout << "SVD solution" << std::endl;
    SVD SVDFilip(AFilip);
    VecDoub aSVDFilip(11);
    SVDFilip.solve(yFilip,aSVDFilip, SVDFilip.eps);
    aSVDFilip.print();

    // Print out machine precision
    std::cout << "\n Machine epsilon (Machine Precision)" << std::endl;
    std::cout << numeric_limits<Doub>::epsilon() << std::endl;

    double test[] = {2.0, 2.0, 1.0};
    double x1Arr[] = {2.0,8.0,4.0,2.0,1.0};
    double x2Arr[] = {1.0,1.0,5.0,7.0,8.0};
    double x3Arr[] = {4.0,-5.0,1.0,-4.0,3.0};
    VecDoub x1(5, x1Arr);
    VecDoub x2(5, x2Arr);
    VecDoub x3(5, x3Arr);

    std::vector<VecDoub> basis = {x1,x2,x3};

    std::cout << "\nCompute orthonormal basis with Gram-Schmidt method" << std::endl;
    std::vector<VecDoub> orth_basis = util::gramSchmidt(basis);

    double bArr[] = {5.0,6.0,1.0,2.0,3.0};
    VecDoub b(5, bArr);

    for (auto &u : orth_basis)
        u.print();

    std::cout << "\nCompute closest point" << std::endl;
    VecDoub b_LS = util::leastSquares(orth_basis, b);
    b_LS.print();

}

/*
 * Examining singular value threshold and impact on coeffiecients
 * Examining error residuals
 * Examining error on the result
 */
void lesson4()
{


    VecDoub xFilip(82); VecDoub yFilip(82);
    util::loadDataset("FilipData.dat",82,xFilip,yFilip);
    MatDoub AFilip = util::getDesignMatrix(xFilip, 11);

    std::cout << "SVD solution" << std::endl;
    SVD SVDFilip(AFilip);
//    VecDoub aSVDFilip(11);
//    SVDFilip.solve(yFilip,aSVDFilip, SVDFilip.eps);
//    aSVDFilip.print();

//    std::cout << "\nSingular values:" << std::endl;
//    SVDFilip.w.print();

//    std::cout << "\nMachine Precision: " << numeric_limits<Doub>::epsilon() <<std::endl;
//    std::cout << "\nReciprocal Condition number: " << SVDFilip.inv_condition() << std::endl;

//    std::cout << "\nResidual error: " << util::residualError(AFilip, aSVDFilip, yFilip) << std::endl;

    std::vector<double> thresholds = {0.01, 0.001, 0.0001, 0.00001, 0.000000001, SVDFilip.eps};
    util::examineResiduals(AFilip, yFilip, 11, thresholds);

    // new design matrix considering errors
    double sigma = 0.00335;
    MatDoub AFilipError = AFilip/sigma;
    //std::cout << "AFilipError" << std::endl;
    //util::print(AFilip);
    //util::print(AFilipError);
    SVD SVDEFilip(AFilipError);

    std::cout << "\nSVD object threshold: " << SVDEFilip.tsh << std::endl;
    std::cout << "\n\n\nError estimates on the solution: " << std::endl;
    VecDoub err_est = util::errorEstimates(AFilipError, 0.00001);
    err_est.print();
    std::cout << "\n\n\nError estimates on the solution with higher threshold: " << std::endl;
    err_est = util::errorEstimates(AFilipError, SVDEFilip.tsh);
    err_est.print();
}

/*
 * Solving non-linear equations
 */
double func(double x)
{
    return x - cos(x);
}




template <class T>
double bisection(T &func, const double x1, const double x2, const double acc)
{

}


void lesson5()
{
    // Solving function
    // x - cos(x) = 0

    double x1 = 0.0;
    double xh = pi/2;
    double acc = std::pow(10,-8);

    // Solve func using Bisection in the interval 0, pi/2

    //std:: cout << rtbis_wp(func,x1,xh,acc) << std::endl;

    //std::cout << zriddr(func, x1, xh, std::pow(10,-16)) << std::endl;


}

void assignment1()
{
    double sigma = 0.1;

    //Load datasets
    VecDoub xd1(500);
    VecDoub yd1(500);
    VecDoub th1d1(500);
    VecDoub th2d1(500);

    util::loadDataset("d1", 500, th1d1, th2d1, xd1, yd1);

    MatDoub Ad1 = as1::getDesignMatrix(th1d1, th2d1);
    VecDoub zd1 = as1::getZ(xd1, yd1);

    SVD SVDd1(Ad1);
    VecDoub qd1(4);
    SVDd1.solve(zd1, qd1, SVDd1.eps);
    std::cout << "d1 estimated parameters" << std::endl;
    qd1.print();
    std::cout << "d1 Singular values: " << std::endl;
    SVDd1.w.print();
    std::cout << "Reciprocal Condition number: " << SVDd1.inv_condition() << std::endl;
    std::cout << "Residual error: " << util::residualError(Ad1, qd1, zd1) << std::endl;

    MatDoub Ad1_err = Ad1/sigma;
    VecDoub err_est_d1 = util::errorEstimates(Ad1_err, SVDd1.tsh);
    err_est_d1.print();

    std::cout << "\n\n";
    //Load datasets
    VecDoub xd2(500);
    VecDoub yd2(500);
    VecDoub th1d2(500);
    VecDoub th2d2(500);

    util::loadDataset("d2", 500, th1d2, th2d2, xd2, yd2);

    MatDoub Ad2 = as1::getDesignMatrix(th1d2, th2d2);
    VecDoub zd2 = as1::getZ(xd2, yd2);

    SVD SVDd2(Ad2);
    VecDoub qd2(4);
    SVDd2.solve(zd2, qd2, SVDd2.eps);
    std::cout << "d2 estimated parameters" << std::endl;
    qd2.print();
    std::cout << "d2 Singular values: " << std::endl;
    SVDd2.w.print();
    std::cout << "Reciprocal Condition number: " << SVDd2.inv_condition() << std::endl;
    std::cout << "Residual error: " << util::residualError(Ad2, qd2, zd2) << std::endl;
    MatDoub Ad2_err = Ad2/sigma;
    VecDoub err_est_d2 = util::errorEstimates(Ad2_err, SVDd2.tsh);
    err_est_d2.print();

    //std::vector<double> thresholds = {0.01, 0.001, 0.0001, 0.00001, 0.000000001, SVDFilip.eps};
    //util::examineResiduals(AFilip, yFilip, 11, thresholds);
}

// Material Constants
const double l7_v = 120; // kg
const double l7_k = 2.5; // meters
const double l7_w = 4.0; // kg/m
const double l7_alpha = 2.0e-7; // kg⁻1

double l7_d = 30; // meters
double l7_n = 5.0;

VecDoub vecfunc(VecDoub_I x)
{
    VecDoub temp(x.size(), 0.0);

    double p = x[0];
    double L = x[1];
    double x_par = x[2];
    double theta = x[3];
    double a = x[4];
    double phi = x[5];
    double L0 = x[6];
    double H = x[7];

    // p
    temp[0] = a * ( cosh( x_par/a ) - 1) - p;
    // L
    temp[1] = 2 * a * sinh( x_par/a ) - L;
    // x (x_par)
    temp[2] = l7_d/2 - l7_k * cos(theta) - x_par;
    // theta
    temp[3] = asin((l7_n-p)/l7_k) - theta;
    // a
    temp[4] = x_par/(asinh(tan(phi))) - a;
    // phi
    temp[5] = atan( tan(theta)/( 1 + l7_v/(l7_w*L0) ) ) - phi;
    // L0
    temp[6] = L/( 1+l7_alpha*H ) - L0;
    // H
    temp[7] = (l7_w*L0) / ( 2*sin(phi) ) - H;

    return temp;
}

//VecDoub vecfunc(VecDoub_I x)
//{
//    VecDoub temp(x.size(), 0.0);

//    double p = x[0];
//    double L = x[1];
//    double x_par = x[2];
//    double theta = x[3];
//    double a = x[4];
//    double phi = x[5];
//    double L0 = x[6];
//    double H = x[7];

//    // p
//    temp[0] = a * ( cosh( x_par/a ) - 1) - p;
//    // L
//    temp[1] = 2 * a * sinh( x_par/a ) - L;
//    // x
//    temp[2] = 30.0/2 - 2.5 * cos(theta) - x_par;
//    // theta
//    temp[3] = asin((5.0-p)/2.5) - theta;
//    // a
//    temp[4] = x_par/(asinh(tan(phi))) - a;
//    // phi
//    temp[5] = atan( tan(theta)/( 1 + 120.0/(4.0*L0) ) ) - phi;
//    // L0
//    temp[6] = L/( 1+(2.0e-7)*H ) - L0;
//    // H
//    temp[7] = (4.0*L0) / ( 2*sin(phi) ) - H;

//    return temp;
//}

void lesson7()
{


    // Start Guesses
    // sg: start guess for that variable
    double sg_phi = pi/6.0; // degrees
    double sg_theta = pi/3.0;
    double sg_x = 13.5; // meters
    double sg_p = 3.0;
    double sg_a = 40.0; //parameter in the catenary equation for the cable
    double sg_H = 100.0; // String tension in the cable
    double sg_L = 30.0; // meters
    double sg_L0 = 30.0;

    double arr_sgx[] = {sg_p, sg_L, sg_x, sg_theta, sg_a, sg_phi, sg_L0, sg_H};

    // Other parameters
    //double d = 30; // meters
    double arr_n[] = {5.0, 2.0, 1.0, 0.5, 0.2, 0.1};

    VecDoub sgx(8,arr_sgx);


    //VecDoub result = vecfunc(sgx);
    //result.print();





    std::cout << "Initial guesses: " << std::endl;
    sgx.print();

    for (auto &sg_n : arr_n)
    {
        l7_n = sg_n;
        std::cout << "\n\n Finding roots for n: " << sg_n << std::endl;
        VecDoub x_res = sgx;
        bool check_var;
        newt(x_res, check_var, vecfunc);
        std::cout << "Estimated parameters: " << std::endl;
        x_res.print();
        std::cout << "Output quantitity check: " << check_var << std::endl;


    }





//    double test[] = {2.0, 2.0, 1.0};
//    double x1Arr[] = {2.0,8.0,4.0,2.0,1.0};
//    double x2Arr[] = {1.0,1.0,5.0,7.0,8.0};
//    double x3Arr[] = {4.0,-5.0,1.0,-4.0,3.0};
//    VecDoub x1(5, x1Arr);
//    VecDoub x2(5, x2Arr);
//    VecDoub x3(5, x3Arr);
}

namespace l8 {

/* Functions whose integrals to be estimated in lesson 8 */
double f1(double x){
    return cos(x*x)*exp(-x);
}

double f2(double x){
    return sqrt(x)*cos(x*x)*exp(-x);
}

double f3(double x){
    return 1.0/sqrt(x)*cos(x*x)*exp(-x);
}

double f4(double x){
    return 1000.0*exp(-1.0/x)*exp(1.0/(1.0-x));
}

double f5(double x){
    return 1.0/sqrt(x)*cos(5*x);
}

}
/* End of namespace l8 */


void lesson8()
{
    std::cout << "First function" << std::endl;
    std::cout << "Rectangular method: " << integral_est::rect(0.0, 1.0, 8, l8::f1) << std::endl;
    std::cout << "Trapezoidal method: " << integral_est::trpz(0.0, 1.0, 8, l8::f1) << std::endl;
    std::cout << "Simpsons method: " << integral_est::simp(0.0, 1.0, 8, l8::f1) << std::endl;

    std::cout << "\n\n";
    std::cout << "Second function" << std::endl;
    std::cout << "Simpsons method: " << integral_est::simp(0.0, 1.0, 8, l8::f2) << std::endl;

    std::cout << "\n\n";
    std::cout << "Third function" << std::endl;
    std::cout << "Rectangle method: " << integral_est::rect(0.0, 1.0, 8, l8::f3) << std::endl;

    std::cout << "\n\n";
    std::cout << "Fourth function" << std::endl;
    std::cout << "Trapezoidal method: " << integral_est::trpz(0.0, 1.0, 8, l8::f4) << std::endl;
    std::cout << "Rectangle method: " << integral_est::rect(0.0, 1.0, 8, l8::f4) << std::endl;
    std::cout << "Simpsons method: " << integral_est::simp(0.0, 1.0, 8, l8::f4) << std::endl;

    util::doSummaryTable(8,2,l8::f1);
    //util::doSummaryTable(2,2,integral_est::rect);
    //util::doSummaryTable(2,2,l8::f1,integral_est::rect<double>);
}

void lesson9()
{
    //std::cout << "First function" << std::endl;
    //double blup = integral_est::de(0.0, 1.0, 4, l8::f5);
    //std::cout << "Double Exponential Rule (DE): " << blup << std::endl;
    //util::doSummaryTable();


    // Previous Mandatory Assignment
    const int T1 = 1000;
    const int T2 = 500;
    const double eps1 = 0.80;
    const double eps2 = 0.60;
    const double sigma = 1.72e-9;
    const double d = 1.00;
    const double w = 1.00;
    const double c11 = eps1*sigma*pow(T1,4);
    const double c12 = (1-eps1);
    const double c21 = eps2*sigma*pow(T2,4);
    const double c22 = (1-eps2);

    // number of iterations
    const int N = 4;

    // right hand side b
    //int * b = new int[]
    VecDoub b(2*N);
    // fill with constants
    for (int i = 0; i < b.size(); i++)
    {
        b[i] = -c11;
        b[i+N] = -c21;
    }

    // unknown parameters x
    VecDoub x(2*N,0.0);

    // step size
    double h = w / (double)(N - 1);

    // x and y
    double x_it = -0.5 * w;
    double y_it = -0.5 * w;
    double scale = 1.0;

    // Design matrix A
    MatDoub A(2*N,2*N,0.0);
    // Determine left top corner of A
    for (int i = 0; i < A.nrows() / 2; i++)
    {
        y_it = -0.5*w;
        x_it += h;

        scale = (i == 0 || i == (A.nrows() / 2 - 1)) ? 0.5 : 1.0;

        for (int j = 0; j < A.ncols() / 2; j++)
        {
            A[i][j] = scale*c12 * h * 0.5 * 1.0 / pow((d*d + (x_it-y_it)*(x_it-y_it)),1.5);
            std::cout << "y_it is: " << y_it << std::endl;
            y_it += h;

        }
    }

    // x and y
    x_it = -0.5 * w;
    y_it = -0.5 * w;
    scale = 1.0;

    // Determine right bottom corner of A
    for (int i = A.nrows() / 2; i < A.nrows(); i++)
    {
        y_it += h;
        x_it = -0.5*w;

        scale = (i == A.nrows() / 2 || i == (A.nrows() - 1)) ? 0.5 : 1.0;

        for (int j = A.ncols() / 2; j < A.ncols(); j++)
        {
            A[i][j] = scale*c22 * h * 0.5 * 1.0 / pow((d*d + (x_it-y_it)*(x_it-y_it)),1.5);
            x_it += h;
        }
    }

    A.print();


}

/*
 * Resources for lesson 10
 */
namespace l10 {

VecDoub f1(double x, VecDoub &y)
{
    VecDoub res(2,0.0);
    res[0]=y[0]*y[1];
    res[1]=-y[0]*y[0];
    return res;
}

VecDoub f2(double x, VecDoub &y)
{
    VecDoub res(2,0.0);
    res[0]=cos(-1 + x + y[0] + 3 * y[1]);
    res[1]=-(y[0]*y[0]) + 2 * sin(y[1]);
    return res;
}

}
/* End of namespace l10 */



/*
 * ODE estimation
 */
void lesson10()
{
    // declaring the function objects to be used for ODE estimation
    ode_est::euler euler_est;
    ode_est::mid mid_est;
    ode_est::leapf leapf_est;
    ode_est::trpz trpz_est;
    ode_est::rk4tho rk4tho_est;

    // testing with at a fixed step-sizel10
    VecDoub y0(2);
    y0[0] = 1.0;
    y0[1] = 1.0;
    est_val l10f1(y0, 0.0, 20.0, 1.0);

//    VecDoub res;
//    double h = 0.01
//    res = euler_est(h, 0, 20, y0, l10::f1);
//    res.print();
//    res = mid_est(h, 0, 20, y0, l10::f1);
//    res.print();
//    res = leapf_est(h, 0, 20, y0, l10::f1);
//    res.print();
//    res = trpz_est(h, 0, 20, y0, l10::f1);
//    res.print();
//    res = rk4tho_est(h, 0, 20, y0, l10::f1);
//    res.print();

    // Make summaries evaluating f1 for 10 iterations of decreasing
    // step-size with the five different ODE estimation methods
    util::summary_table_ode_est(10, l10f1, l10::f1, euler_est);
    util::summary_table_ode_est(10, l10f1, l10::f1, mid_est);
    util::summary_table_ode_est(10, l10f1, l10::f1, leapf_est);
    util::summary_table_ode_est(10, l10f1, l10::f1, trpz_est);
    util::summary_table_ode_est(10, l10f1, l10::f1, rk4tho_est);


}

void lesson11()
{
    ode_est::mid mid_est;

    VecDoub y0(2);
    y0[0] = 1;
    y0[1] = 0;
    est_val l11f2(y0, 0.0, 1.0, 1.0);

    util::summary_table_ode_est(10, l11f2, l10::f2, mid_est);
}

#endif // LESSONS_H
