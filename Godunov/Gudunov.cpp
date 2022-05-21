/**
* @mainpage  Riemann Solvers
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
 * @file Gudunov.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief The second version of Godunov scheme in which the CFL like condition is dt \leq dx/S_{max}.
 * The scheme can be written in a conservative form, where the numerical flux f_{i+1/2} is the exact flux 
 * evaluated with the exact Riemann solver.
 * @version 1.01
 * @date 2022-05-21
 * 
 * @copyright Copyright (c) 2022 R.Y. Han, All Rights Reserved.
 * 
 */

#include <iostream>
#include <cmath>
#include <fstream>
#include "configration.h"
#include "sampling.h"
using namespace std;

/**
 * @brief The exact flux in Euler equations with conservative variables;
 * 
 * @param states : (\rho, u, p)
 * @return f(\vector{u})
 * @note The input states are the primary variables (rho, u, p), but the flux is
 * evaluated with the conservative variables (rho, rho*u, E).
 */
double* flux(double *states);

/**
 * @brief Initialize the primary states.
 * 
 * @param un : The vector to storage numerical results.
 * @note Only initialize the PRIMARY variables.
 */
void Initial(double **un);

/**
 * @brief Godunov flux
 * 
 * @param un : Numerical results at step t_n
 * @param j : The index of the local Riemann problem
 * @return f_{j+1/2}
 * @note The input un are the primary variables (rho, u, p) but the flux is evaluated with conservative variables
 */
double* Godunovflux(double *states, int j);

/**
 * @brief Choose the time step to ensure CFL like condition.
 * 
 * @param states : Primary variables.
 * @return time step.
 * @note Input states are primary variables.
 */
double choosedt(double **states);

/**
 * @brief The conservative schemes, in which the numerical flux is Godunov flux.
 * 
 * @return Numerical reusults at step t_{n+1}.
 */
double** conservativeScheme();

/**
 * @brief Plot the numerical results at the final time T.
 * 
 * @param states : Numerical results of primary variables.
 */
void plot(double **states);

int main()
{
    int j;
    double **states = new double*[N];
    for (j=0; j<N; j++)
    {
        states[j] = new double [3];
    }
    states = conservativeScheme();

    plot(states);

    system("pause");
}

double* flux(double *states)
{
    double *ans = new double [3];
    double a1, a2, a3;
    a1 = states[0] * states[1]; /** rho * u */

    a2 = states[0] * states[1] * states[1] + states[2]; /** rho * u^2 + p */

    a3 = states[2] / ((gamma-1)*states[0]); /** e */
    a3 = states[0] * (0.5 * states[1]*states[1] + a3);  /** E */
    a3 = states[1] * (a3 + states[2]);  /** u(E+p) */

    ans[0] = a1;
    ans[1] = a2;
    ans[2] = a3;
    return ans;
}

void Initial(double **un)
{
    double x;
    int j;
    for (j=0; j<N; j++)
    {
        x = a + (j+0.5) * h;    /** x is the central point of cell j */
        if (x <= InitialDiscontinuity)
        {
            un[j][0] = rhol_ini;
            un[j][1] = ul_ini;
            un[j][2] = pl_ini;
        }
        else{
            un[j][0] = rhor_ini;
            un[j][1] = ur_ini;
            un[j][2] = pr_ini;
        }
    }
}

double* Godunovflux(double **states, int j)
{
    double *statesl, *statesr;
    double* ans;
    if (j == -1)
    {
        statesl = states[0];
        statesr = states[0];
    }
    else if (j == N-1)
    {
        statesl = states[N-1];
        statesr = states[N-1];
    }
    else{
        statesl = states[j];
        statesr = states[j+1];
    }
    
    ans = sampling(0,T,statesl,statesr);  /** u_{i+1/2}(0) of primary variables */
    ans = flux(ans);

    return ans;
}

double choosedt(double **states)
{
    double dt, smax, temp;

    smax = 0;
    for (int j=0; j<N; j++)
    {
        temp = sqrt((gamma-1)* states[j][2] / states[j][0]);    /** a_j */
        temp = abs(states[j][1]) + temp;    /** |u_j| + a_j */

        if (temp > smax)
        {
            smax = temp;
        }
    }

    dt = CFL * h / smax;
    return dt;
}

double** conservativeScheme()
{
    double **states = new double* [N];   /** primary variables */
    double **un = new double* [N];   /** conservative variables */
    double **temp = new double* [N];
    double *f1 = new double [3];
    double *f2 = new double [3];
    int i,j;
    double t = 0, dt;
    for (j=0; j<N; j++)
    {
        states[j] = new double [3];
        un[j] = new double [3];
        temp[j] = new double [3];
    }

    /** First initialize the variables (primary and conservative) */
    Initial(states);
    for (j=0; j<N; j++)
    {
        un[j][0] = states[j][0];    /** rho */
        un[j][1] = states[j][0] * states[j][1]; /** rho * u */
        un[j][2] = states[j][2] / ((gamma-1) * states[j][0]);  /** e */
        un[j][2] = states[j][0] * ( 0.5 * states[j][1] * states[j][1] + un[j][2]);  /** E */
    }

    /** solve the equation */
    while(t<T)
    {
        /** choose dt */
        if ( t + dt >T)
        {
            dt = T-t;
        }
        else{
            dt = choosedt(states);
        }
        t = t + dt;
        cout<<t<<"\t"<<dt<<endl;

        /** time evolution */
        temp = un;

        for (j=0; j<N; j++)
        {
            f1 = Godunovflux(states,j-1);
            f2 = Godunovflux(states,j);
            for (i=0; i<3; i++)
            {
                un[j][i] = temp[j][i] + dt * ( f1[i] - f2[i] )/h;
            }
        }

        for (j=0; j<N; j++)
        {
            states[j][0] = un[j][0];  /** rho */
            states[j][1] = un[j][1] / un[j][0]; /** u */
            states[j][2] = un[j][2] / un[j][0]- 0.5*un[j][1] * un[j][1];   /** e */
            states[j][2] = states[j][2] * (gamma-1) * un[j][0]; /** p */
        }

    }

    for (j=0; j<N; j++)
    {
        delete[] un[j];
        delete[] temp[j];
    }
    delete[] un, temp;

    return states;
}

void plot(double **states)
{
    double et1, et2, *temp, *sl, *sr, *x;
    temp = new double[3];
    sl = new double[3];
    sr = new double[3];
    x = new double [N];
    int j;

    sl[0] = rhol_ini;
    sl[1] = ul_ini;
    sl[2] = pl_ini;
    sr[0] = rhor_ini;
    sr[1] = ur_ini;
    sr[2] = pr_ini;

    for (j=0; j<N; j++)
    {
        x[j] = a + (j+0.5) * h;
    }

    const char* fnrho = "F:\\C++code\\RiemannSolver\\Godunov\\results\\density.plt";
    const char* fnu = "F:\\C++code\\RiemannSolver\\Godunov\\results\\velocity.plt";
    const char* fnp = "F:\\C++code\\RiemannSolver\\Godunov\\results\\pressure.plt";
    const char* fne = "F:\\C++code\\RiemannSolver\\Godunov\\results\\internal_energy.plt";

    remove(fnrho);
    remove(fnu);
    remove(fnp);
    remove(fne);

    fstream fr, fu, fp, fe;

    fr.open(fnrho, ios :: out | ios :: app);
    fu.open(fnu, ios :: out | ios :: app);
    fp.open(fnp, ios :: out | ios :: app);
    fe.open(fne, ios :: out | ios :: app);
    fr<<"VARIABLES = X, rho_h, rho"<<endl;
    fu<<"VARIABLES = X, u_h, u"<<endl;
    fp<<"VARIABLES = X, p_h, p"<<endl;
    fe<<"VARIABLES = X, e_h, e"<<endl;
    for (j=0; j<N; j++)
    {
        temp = sampling(x[j]-InitialDiscontinuity,T,sl,sr);
        fr<<"\t"<<x[j]<<"\t"<<states[j][0]<<"\t"<<temp[0]<<endl;
        fu<<"\t"<<x[j]<<"\t"<<states[j][1]<<"\t"<<temp[1]<<endl;
        fp<<"\t"<<x[j]<<"\t"<<states[j][2]<<"\t"<<temp[2]<<endl;

        et1 = states[j][2] / ((gamma - 1)* states[j][0]);
        et2 = temp[2] / ((gamma-1)*temp[0]);
        fe<<"\t"<<x[j]<<"\t"<<et1<<"\t"<<et2<<endl;
    }
    fr.close();
    fu.close();
    fp.close();
    fe.close();

    delete[] sl, sr, temp;
}