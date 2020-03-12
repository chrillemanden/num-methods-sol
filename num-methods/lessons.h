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

    std::cout << zriddr(func, x1, xh, std::pow(10,-16)) << std::endl;


}

void assignment1()
{
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

    //std::vector<double> thresholds = {0.01, 0.001, 0.0001, 0.00001, 0.000000001, SVDFilip.eps};
    //util::examineResiduals(AFilip, yFilip, 11, thresholds);


}

#endif // LESSONS_H
