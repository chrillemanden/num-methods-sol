#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <iostream>
#include <string>
#include <memory>



/* Namespace containing:
 * Function objects for approximating integrals
 *  - Extended Midpoint (Rectangle) Method
 *  - Trapezoidal Method
 *  - Simpsons Method
 *
 * All methods follow this structure:
 * @param lb: lower bound for the integral
 * @param ub: upper bound for the integral
 * @param h: step size when running the algorithm
 * @param func: function to estimate integral for
 * @return estimated integral value
 */
namespace int_est {

struct InteStartCond
{
    double lb; // lower bound
    double ub; // upper bound
    double h; // initial step size
    InteStartCond(double lb, double ub, double h) : lb(lb), ub(ub), h(h) {}
};

/*
 * The function object definitions are derived
 * from this definition
 */
struct IntegralEst {
    std::string const desc;
    int order;
    int f_comp;
    virtual double operator()(double lb, double ub, double h, double (*func) (double)) = 0;
    virtual ~IntegralEst() noexcept = default;
protected:
    IntegralEst(std::string desc, int order, int f_comp) noexcept
        : desc(move(desc)), order(order), f_comp(f_comp) {}
};

struct Rect : public IntegralEst {
    Rect() : IntegralEst("Extended Midpoint (Rectangle) Method", 2, 0) {}
    double operator()(double lb, double ub, double h, double (*func) (double)) override
    {
        f_comp = 0;

        // Initialise estimate
        double Ah = 0.0;

        // Accumulate estimate
        for (double mid = lb + 0.5 * h; mid < ub; mid+=h)
        {
            Ah += h*func(mid);
            f_comp++;
        }
        return Ah;
    }
};

struct Trpz : public IntegralEst {
    Trpz() : IntegralEst("Trapezoidal Method", 2, 0) {}
    double operator() (double lb, double ub, double h, double (*func) (double)) override
    {
        f_comp = 1;

        // Initialise the estimate
        double Ah = 0.5 * func(lb) * h;

        // Accumulate the estimate
        for (double mid = lb + h; mid < ub; mid+=h)
        {
            Ah += h * func(mid);
            f_comp++;
        }

        Ah += 0.5 * func(ub) * h;
        f_comp++;

        return Ah;
    }
};

struct Simp : public IntegralEst {
    Simp() : IntegralEst("Simpsons Method", 4, 1) {}
    double operator() (double lb, double ub, double h, double (*func) (double)) override
    {
        f_comp = 1;
        // Initialise the estimate
        double Ah = 1.0/3.0 * func(lb) * h;

        // Accumulate the estimate
        int i = 1;
        for (double mid = lb + h; mid < ub; mid+=h)
        {
            Ah += (2 + (i % 2) * 2)/3.0 * h * func(mid);
            i++;
            f_comp++;
        }

        // Evaluate at upper bound
        Ah += 1.0/3.0 * func(ub) * h;
        f_comp++;

        return Ah;
    }
};

template <typename P>
double de(double lb, double ub, int N, P &func)
{
    // infinity approximation (5.0 is very close to infinity!)
    double T = 5.0;

    // step_size
    double h = (T - (-T)) / N;
    std::cout << "h is " << h << std::endl;

    // initial_estimate
    double Ah = 0.0;
    for(double t = -T; t <= T; t+=h)
    {
        std::cout << "t: " << t << std::endl;
        double x = 0.5 * (ub + lb) + 0.5 * (ub - lb)*tanh(sinh(t));

        std::cout << "x: " << x << std::endl;
        double dxdt = exp(-1.0*exp(abs(t)));
        std::cout << "dxdt: " << dxdt << std::endl;
        double q = exp(-2 * sinh(t));
        std::cout << "q: " << q << std::endl;
        double dxdt2 = 2.0*(ub-lb)*q/((1+q)*(1+q))*cosh(t);
        std::cout << "dxdt2: " << dxdt2 << std::endl;
        Ah += func(x) * dxdt2;
    }
    return Ah * h;
}

}
/* End of namespace integral_est */

#endif // INTEGRAL_H
