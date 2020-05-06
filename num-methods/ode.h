#ifndef ODE_H
#define ODE_H

#include <iostream>
#include <string>
//#include <memory>

// Numerical Recipes inlcludes
#include "nr3.h"
#include "tridag.h"

#include "opol.h"

/*
 * Object for trapezoidal method for initial value problem estimation
 */
struct Phi_s
{
    VecDoub yn;
    double xn;
    VecDoub f;
    double h;
    VecDoub (*func)(double, VecDoub&);
    //Phi_s(VecDoub yn, VecDoub f, double h, VecDoub &func) : yn(yn), f(f), h(h), func(func) {}
    VecDoub operator() (VecDoub y)
    {
        return y - yn - (f + func(xn + h, y)) * 0.5 * h;
        //return y - yn - f * 0.5 * h;
        //return y;
    }
    void setParams(double xnn, VecDoub &ynn, VecDoub &fn)
    {
        xn = xnn;
        yn = ynn;
        f = fn;
    }
};


// Bounded Value Problem Starting Conditions
//template <typename T1, typename T2, typename T3>
struct BVPsc
{
    double (*F) (double, double, double);
    double (*Fy) (double, double, double);
    double (*Fyp) (double, double, double);
//    T1 &F;
//    T2 &Fy;
//    T3 &Fyp;
    double alpha;
    double beta;
    double a;
    double b;
    double trsh; // threshold for the newton estimation
    double trg_x; // the x-value where we want to find y
    BVPsc(double (*F) (double, double, double), double (*Fy) (double, double, double), double (*Fyp) (double, double, double), double alpha, double beta, double a, double b, double trsh, double trg_x)
    //rel_s(T1 &F, T2 &Fy, T3 Fyp, double alpha, double beta, double a, double b)
        : F(F), Fy(Fy), Fyp(Fyp), alpha(alpha), beta(beta), a(a), b(b), trsh(trsh), trg_x(trg_x) {}
};

struct JacPhi
{
    VecDoub ja; // lower diagonal
    VecDoub jb; // diagonal
    VecDoub jc; // upper diagonal
    VecDoub phi; // right hand side
    JacPhi(int N)
        : ja(VecDoub(N-1, 0.0)), jb(VecDoub(N-1, 0.0)), jc(VecDoub(N-1, 0.0)), phi(VecDoub(N-1, 0.0)) {}
};


namespace ode_est {

struct euler
{
    std::string desc = "1st Order Runge-Kutta (Euler) Method";
    int order = 1;
    int f_comp = 0;
    VecDoub operator() (double h, double x0, double x, const VecDoub &y0, VecDoub (*func)(double, VecDoub&))
    {
        f_comp = 0;
        VecDoub yn = y0;

        for (double xn = x0; xn < x; xn += h)
        {
            yn = yn + func(xn, yn) * h;
            f_comp++;
        }
        return yn;
    }
};



struct mid
{
    std::string desc = "2nd Order Runge-Kutta (Midpoint) Method";
    int order = 2;
    int f_comp = 0;
   // mid() {}
    VecDoub operator() (double h, double x0, double x, const VecDoub &y0, VecDoub (*func)(double, VecDoub&))
    {
        f_comp = 0;
        VecDoub yn = y0;
        VecDoub k1;
        VecDoub k2;
        VecDoub tmp;

        for (double xn = x0; xn < x; xn += h)
        {
            k1 = func(xn,yn)*h;
            tmp = yn + k1 * 0.5;
            k2 = func(xn + 0.5 * h,tmp)*h;
            yn = yn + k2;

            f_comp+=2;
        }
        return yn;
    }
};



struct trpz
{
    std::string desc = "Trapezoidal Method";
    int order = 2;
    int f_comp = 0;
   // trpz() {}

    VecDoub operator() (double h, double x0, double x, const VecDoub &y0, VecDoub (*func)(double, VecDoub&))
    {
        f_comp = 0;
        VecDoub yn = y0;
        VecDoub f1;
        VecDoub y_euler;
        bool check_var;

        //Phi_s phi(yn,f1,h,func);
        Phi_s phi;
        phi.h = h;
        phi.func = func;

        for (double xn = x0; xn < x; xn += h)
        {
            f1 = func(xn, yn);
            y_euler = yn + f1 * h;
            phi.setParams(xn, yn, f1);
            yn = y_euler;
            newt(yn, check_var, phi);
            f_comp += 3;
        }
        return yn;
    }
};


struct rk4tho
{
    std::string desc = "4th Order Runge-Kutta Method";
    int order = 4;
    int f_comp = 0;
   // rk4tho() {}

    VecDoub operator() (double h, double x0, double x, const VecDoub &y0, VecDoub (*func)(double, VecDoub&))
    {
        f_comp = 0;
        VecDoub yn = y0;
        VecDoub k1;
        VecDoub k2;
        VecDoub k3;
        VecDoub k4;
        VecDoub tmp;


        for (double xn = x0; xn < x; xn += h)
        {
            k1 = func(xn, yn)*h;
            tmp = yn + k1 * 0.5;
            k2 = func(xn + 0.5 * h, tmp)*h;
            tmp = yn + k2 * 0.5;
            k3 = func(xn + 0.5 * h, tmp)*h;
            tmp = k3 +yn;
            k4 = func(xn + h, tmp)*h;
            yn = yn + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

            f_comp += 4;
        }
        return yn;
    }

};



struct leapf
{
    std::string desc = "Two-step Leap-frog Method";
    int order = 2;
    int f_comp = 1;
   // leapf() {}
    VecDoub operator() (double h, double x0, double x, const VecDoub &y0, VecDoub (*func)(double, VecDoub&))
    {
        f_comp = 1;
        VecDoub ynpre = y0;
        VecDoub yn = ynpre + func(x0, ynpre)*h;
        VecDoub yntmp;

        for (double xn = x0+h; xn < x; xn += h)
        {
            yntmp = yn;
            yn = ynpre + func(xn, yn) * 2 * h;
            ynpre = yntmp;
            f_comp++;
        }
        return yn;
    }
};


struct BoundValueEst {
    std::string const desc;
    int order;
    int f_comp;
    virtual double operator()(BVPsc &p, int N) = 0;
    virtual ~BoundValueEst() noexcept = default;
protected:
    BoundValueEst(std::string desc, int order, int f_comp) noexcept
        : desc(move(desc)), order(order), f_comp(f_comp) {}
};


void calc_tri_jacobian_and_phi(BVPsc &p, JacPhi &J, const std::vector<double> &x, const std::vector<double> &y, int N)
{
    //step size
    double h = (p.b-p.a)/(double)N;

    // J(1,1) & J(1,2)
    double yp = (y[1]-p.alpha) / (2*h);
    J.jb[0] = 2.0+h*h*p.Fy(yp,y[0],x[0]);
    J.jc[0] = -1.0+0.5*h*p.Fyp(yp,y[0],x[0]);

    // phi(1)
    J.phi[0] = -2.0*y[0]+y[1]-h*h*p.F(yp,y[0],x[0]) + p.alpha;

    for (int i = 1; i < N - 2; i++)
    {
        // J(i,i-1) & J(i,i) & J
        yp = (y[i+1] - y[i-1]) / (2.0*h);
        J.ja[i] = -1.0-0.5*h*p.Fyp(yp,y[i],x[i]);
        J.jb[i] = 2.0 + h*h*p.Fy(yp,y[i],x[i]);
        J.jc[i] = -1.0+0.5*h*p.Fyp(yp,y[i],x[i]);

        // phi(i)
        J.phi[i] = y[i-1] - 2.0*y[i] + y[i+1] - h*h*p.F(yp,y[i],x[i]);
    }

    // J(N-1,N-2) & J(N-1,N-1)
    yp = (p.beta-y[N-3]) / (2.0*h);
    J.ja[N-2] = -1.0-0.5*h*p.Fyp(yp,y[N-2],x[N-2]);
    J.jb[N-2] = 2.0+h*h*p.Fy(yp,y[N-2],x[N-2]);

    // phi(N-1)
    J.phi[N-2] = y[N-3] - 2*y[N-2] - h*h*p.F(yp,y[N-2],x[N-2]) + p.beta;
}

struct FinDiff : public BoundValueEst {
    FinDiff() : BoundValueEst("Finite Difference (Relaxation) Method", 2, 0) {}
    double operator ()(BVPsc &p, int N) override
    {
        int max_iter = 100;
        f_comp = 0;

        //step size
        double h = (p.b-p.a)/(double)N;

        // index of y(est_x)
        int xi = p.trg_x/(p.b-p.a)*N - 1;
        //std::cout << "Index is: " << xi << std::endl;

        // xi and yi
        std::vector<double> x;
        std::vector<double> y;
        int i;
        for (i = 1; i < N; i++)
        {
            x.push_back(p.a+i*h);
            y.push_back(p.alpha + (double)i/(double)N*(p.beta-p.alpha));
        }

        // Newton Estimation
        i = 1;
        double dyxi = 1.0;
        while(abs(dyxi) > p.trsh && i++ <= max_iter)
        {
            // Compute step
            JacPhi J(N);
            calc_tri_jacobian_and_phi(p, J, x, y, N);
            f_comp += (N-1)+(N-3)*3+4;

            // Solve step
            VecDoub dy(N-1, 0.0);
            tridag(J.ja, J.jb, J.jc, J.phi, dy);

            // Update step
            y = y + dy;

            dyxi = dy[xi];
        }
        return y[xi];
    }
};


}
/* End of namespace ode_est */

#endif // ODE_H
