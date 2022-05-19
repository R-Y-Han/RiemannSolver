/**
* @mainpage  Exact Riemann Solver
* <table>
* <tr><th>Project  <td> Riemann Solvers and Numerical Methods for Fluid Dynamics
* <tr><th>Author   <td> R.Y. Han
* <tr><th>Institution   <td> IAPCM
* </table>
* @section   Project details
* Code the main Riemann solvers introduced in the book
* Riemann Solvers and Numerical Methods for Fluid Dynamics from E.F. Toro. 
* Mainly for 1-dimension Euler equations.
* 
**********************************************************************************
*/

/**
 * @file UpWindScheme.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief Four classical upwind schemes for LINEAR hyperbolic equations
 * @version 1.01
 * @date 2022-05-17
 * 
 * @details Used four classical 1st order upwind schemes to solve the linear 
 * hyperbolic equation: u_t + a u_x = 0, in which a is a constant. 
 * The four schemes are: Godunov scheme, Lax-Friedrichs scheme,
 * Lax-Wendroff scheme and Warming-Beam scheme.         \n 
 * @htmlonly 
 * <span style="font-weight: bold">History</span> 
 * @endhtmlonly 
 * Project|Version|Auther|Date|Describe
 * -------------|------|----|------|-------- 
 * Riemann Solvers \n and Numerical Methods \n  for Fluid Dynamics|V1.01|R.Y. Han|5.17.2022|4 upwind schemes for linear hyperbolic equations
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <iostream>
#include <cmath>
#include <fstream>
#include "configration.h"
using namespace std;

/**
 * @name Gauss-Lobatto quadrature
 * @brief Gauss-Lobatto quadrature points and coefficients. Used to compute the cell average.
 * @{
 */
const double lobattopoint[5] = {-1.0, -0.6546536707079771, 0, 0.6546536707079771, 1.0};
const double lobattoco[5] = {0.1, 0.5444444444444444, 0.7111111111111111, 0.5444444444444444, 0.1};
/** @} */

/**
 * @brief Godunov flux (at the point x_{i+1/2})
 * 
 * @param ul : value at the cell i
 * @param ur : value at the cell i+1
 * @return f_{i+1/2}^{god}
 * @note Only for linear equations
 */
double Godunov(double ul, double ur);

/**
 * @brief Lax-Fridrichs flux( at the point x_{i+1/2})
 * 
 * @param ul : value at the cell i
 * @param ur : value at the cell i+1
 * @return f_{i+1/2}^{fl}
 */
double LFflux(double ul, double ur);

/**
 * @brief Lax-Wendroff flux ( at the point x_{i+1/2})
 * 
 * @param ul : value at the cell i
 * @param ur : value at the cell i+1
 * @return f_{i+1/2}^{LW}
 */
double LWflux(double ul, double ur);

/**
 * @brief Warming and Beam flux (at the point x_{i+1/2})
 * 
 * @param ul : value at the cell i
 * @param ur : value at the cell i+1
 * @return f_{i+1/2}^{WB}
 */
double WBflux(double ul, double ur);

/**
 * @brief The general conservative scheme
 * in which the numerical flux function only have two parameters: 
 * f_{i+1/2} (u_{i-1}, u_i)
 * 
 * @param flux The numerical flux function, can be chosen as godunov flux or other censervative fluxes
 * @return the numerical results at the final time T 
 */
double* conservativeScheme(double (*flux)(double, double));

/**
 * @brief Plot the results of the numerical solution and the exact solution
 * 
 * @param uh : the numerical solution solved by different schemes
 * @param fn : filename
 */
void plot(double* uh, const char *fn);


int main()
{
    double* god = new double [n];
    double* lf = new double [n];
    double* lw = new double [n];
    double* wb = new double [n];

    const char* fngod = "F:\\C++code\\RiemannSolver\\UpWindScheme\\results\\Godunov.plt";
    const char* fnlf = "F:\\C++code\\RiemannSolver\\UpWindScheme\\results\\Lax-Fridrichs.plt";
    const char* fnlw = "F:\\C++code\\RiemannSolver\\UpWindScheme\\results\\Lax-Wendroff.plt";
    const char* fnwb = "F:\\C++code\\RiemannSolver\\UpWindScheme\\results\\Warming-Beam.plt";
    god = conservativeScheme(Godunov);
    lf = conservativeScheme(LFflux);
    lw = conservativeScheme(LWflux);
    wb = conservativeScheme(WBflux);
    plot(god,fngod);
    plot(lf,fnlf);
    plot(lw,fnlw);
    plot(wb,fnwb);
    system("pause");
}


double Godunov(double ul, double ur)
{
    double ans;

    if (a>=0)
    {
        ans = f(ul);
    }
    else
    {
        ans = f(ur);
    }
    return ans;
}

double LFflux(double ul, double ur)
{
    double ans;
    ans = (1+c)/(2*c) * f(ul) + (c-1)/(2*c) * f(ur);
    return ans;
}

double LWflux(double ul, double ur)
{
    double ans;
    ans = (1+c)* ul /2 + (1-c) * ur /2;
    ans = f(ans);
    return ans;
}

double WBflux(double ul, double ur)
{
    double ans;
    ans = (c-1)*f(ul)/2 + (3-c) * f(ur) /2;
    return ans;
}

double* conservativeScheme(double (*flux)(double, double))
{
    double* uh = new double [n];
    double* temp = new double [n];
    int i, j, conter=0;
    double xi, ul, ur;

    //Firstly initialize uh
    for (i=0; i<n; i++)
    {
        uh[i] = 0;
        for (j=0; j<5; j++)
        {
            xi = h * (lobattopoint[j] + 1) / 2.0 + (i * h + left);
            uh[i] = uh[i] + u_0(xi) * lobattoco[j];
        }
        uh[i] = uh[i] / 2.0;
    }

    while(conter < nt)
    {
        conter++;
        
        //restore the results of last time step in temp
        for (i=0; i<n; i++)
        {
            temp[i] = uh[i];
        }

        if (i == 0)
        {
            ul = uh[0];
        }
        else{
            ul = uh[i-1];
        }
        if (i == n-1)
        {
            ur = uh[n-1];
        }
        else{
            ur = uh[i+1];
        }
        
        uh[i] = temp[i] + (dt/h) * (flux(ul, temp[i]) - flux(temp[i],ur));
    }

    delete[] temp;
    return uh;
}

void plot(double* uh, const char* fn)
{
    int i;
    double temp;
    fstream ff;

    remove(fn);

    ff.open(fn, ios :: out | ios :: app);
    ff<<"VARIABLES="<<"X"<<","<<"uh"<<","<<"u"<<endl;
    for (i=0; i<n; i++)
    {
        temp = (i+0.5)*h + left + a * T;
        ff<<"\t"<<temp<<"\t"<<uh[i]<<"\t"<<u_exact(temp,T)<<endl;
    }
    ff.close();
}