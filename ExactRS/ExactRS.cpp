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
* @file [ExactRS.cpp]
* @brief the exact Riemann Solver for the 1-dimension Euler equations.
* @author R.Y. Han,  IAPCM
* @email: hanruoyu21@gscaep.ac.cn
* @date 5.14.2022
* @version v1.01
*
* @details To compute the states in the star region, and identify two nonlinear waves.   \n 
* Plot the result at a given finish time T, in a region I = [a,b].             \n 
* @htmlonly 
* <span style="font-weight: bold">History</span> 
* @endhtmlonly 
* Project|Version|Auther|Date|Describe
* -------------|------|----|------|-------- 
* Riemann Solvers \n and Numerical Methods \n  for Fluid Dynamics|V1.01|R.Y. Han|5.14.2022|sampling with Exact Riemann Solver
* @date 2022-05-14
* 
* @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
*/  

#include <iostream>
#include <cmath>
#include <fstream>
#include "configuration.h"
using namespace std;



/**
 * @brief The function to compute the relation of pressure across the left nonlinear wave.
 * 
 * @param[double] p : Pressure input.  
 * @return double : \n 
 * The input pressure p is larger than the pressure pl in the left region,
 * then use the function of a shock, else use the function of a rarefaction. Return the computed value.
 */
double f_l(double p);

/**
 * @brief The function to compute the relation of pressure across the right nonlinear wave.
 * 
 * @param[double] p : Pressure input
 * @return double : \n 
 * If the input pressure p is larger than the pressure pl in the right region,
 * then use the function of a shock, else use the function of a rarefaction. Return the computed value.
 */
double f_r(double p);

/**
 * @brief f(p) = f_l(p) + f_r(p) + [u]
 * 
 * @param[double] p : Pressure input.
 * @return double : \n 
 * When p = p*, f(p*) = 0.
 */
double f(double p);

/**
 * @brief The derivative function of f_l
 * 
 * @param[doule] p : Pressure input 
 * @return double 
 */
double f_ldp(double p);

/**
 * @brief The derivative function of f_r
 * 
 * @param[double] p : Pressure input 
 * @return double 
 */
double f_rdp(double p);

/**
 * @brief Compute the pressure in the star region using Newton's iterative methods.
 * p* satisfies f(p*) = 0.
 * 
 * @return [double] p* : The pressure in the star region.
 */
double getpstar();

/**
 * @brief Compute the velocity in the star region.
 * 
 * @param[double] pstar : the pressure in the star region 
 * @return[double] u* : the velocity in the star region
 * @warning During iteration the parameter temp may be negative, choose the suitable initial value
 * of TRRS or mean value to avoid error.
 */
double getustar(double pstar);

/**
 * @brief Identify the type of two nonlinear waves.
 * @details 
 * -# If p* > pl, the left wave is a shock;
 * -# If p* < pl, the left wave is a rarefaction;
 * -# If p* > pr, the right wave is a shock;
 * -# If p* < pr, the right wave is a rarefaction.
 * 
 * @param[double] p : pressure input, p = pl or pr
 * @return wavetype : 
 * -# When p = pl, return the type of left wave;
 * -# When p = pr, return the type of right wave.
 */
wavetype identifyNonlinearWave(double p);

/**
 * @brief Compute the density between left wave and the contact discontinuity.
 * 
 * @param pstar : the pressure in the star region
 * @param leftwave : the type of the left wave
 * @return double : the density between the left wave and the contact discontinuity.
 */
double getrhostarl(double pstar, wavetype leftwave);

/**
 * @brief Compute the density between right wave and the contact discontinuity.
 * 
 * @param pstar  : the pressure in the star region
 * @param rightwave : the type of the right wave
 * @return double : the density between the right wave and the contact discontinuity.
 */
double getrhostarr(double pstar, wavetype rightwave);

/**
 * @brief Plot the picture of the exact solutions to the Euler equations.
 * 
 * @param leftwave : type of the left wave
 * @param rightwave : type of the right wave
 * @param dll : the distance travelled of the head of the left wave
 * @param dlr : the distance travelled of the tail of the left wave
 * @param dmid : the distance travelled of the contact discontinuity
 * @param drl : the distance travelled of the tail of the right wave
 * @param drr : the distance travelled of the head of the right wave
 * @param rhostarl : density between left region and contact discontinuity
 * @param rhostarr : density between contact discontinuity and right region
 * @param ustar : velocity in the star region
 * @param pstar : pressure in the star region
 */
void plot(wavetype leftwave, wavetype rightwave, double dll, double dlr, double dmid, double drl, double drr, double rhostarl, double rhostarr, double ustar, double pstar);

/**
 * @brief Solve the Euler equations exactly. \n 
 * Compute the density, velocity, pressure and internal energy in the star region,
 * the distance travelled of the three waves.
 * 
 */
void exactRS();

/**
 * @brief Get the required states at the given location x at the final time T.
 * 
 * @param x : given location
 * @return double* : a vector of 4 components: density, velocity, pressure, internal energy.
 */
double* sampling(double x);


int main()
{
    exactRS();
    system("pause");
}


double f_l(double p)
{
    double ans;
    if (p > pl)
    {
        ans = Al / (p + Bl);
        ans = sqrt(ans);
        ans = ans * (p - pl);
    }
    else{
        double temp;
        temp = (gamma - 1)/(2 * gamma);

        ans = p / pl;
        ans = pow(ans, temp);
        ans = 2 * al * (ans - 1)/(gamma - 1);
    }
    return ans;
}

double f_r(double p)
{
    double ans;
    if ( p > pr)
    {
        ans = Ar / ( p + Br );
        ans = sqrt(ans);
        ans = ans * (p - pr);
    }
    else{
        double temp;
        temp = (gamma - 1) / (2 * gamma);

        ans = p / pr;
        ans = pow(ans, temp);
        ans = 2 * ar * (ans - 1)/ (gamma - 1);
    }
    return ans;
}

double f(double p)
{
    double ans;
    ans = f_l(p) + f_r(p) + (ur - ul);
    return ans;
}

double f_ldp(double p)
{
    double ans;

    if (p > pl)
    {
        ans = 1 - (p - pl)/(2 * (Bl + p));
        ans = ans * sqrt(Al/(Bl + p));
    }
    else{
        double temp;
        temp = (gamma + 1)/(2 * gamma);

        ans = p/pl;
        ans = pow(ans,-temp);
        ans = ans / (rhol * al);
    }
    return ans;
}

double f_rdp(double p)
{
    double ans;

    if (p > pr)
    {
        ans = 1 - (p - pr)/(2 * (Br + p));
        ans = ans * sqrt(Ar/(Br + p));
    }
    else{
        double temp;
        temp = (gamma + 1)/(2 * gamma);

        ans = p/pr;
        ans = pow(ans,-temp);
        ans = ans / (rhor * ar);
    }
    return ans;
}

double getpstar()
{
    double pstar, temp;
    double fp, fdp, e, time;
    e = 1e-6;
    time = 1;
    //pstar = al+ar - 0.5*(gamma-1)*(ur-ul);
    //pstar = pstar / (pow(al/pl,(gamma-1)/(2*gamma)) + pow(ar/pr,(gamma-1)/(2*gamma)));
    //pstar = pow(pstar,2*gamma/(gamma-1));
    pstar = 0.5 * (pl + pr);
    temp = 0;
    while(abs(pstar - temp)>e)
    {
        time++;
        if (time > 1e6)
        {
            break;
        }

        temp = pstar;
        fp = f(temp);
        fdp = f_ldp(temp) + f_rdp(temp);
        pstar = temp - fp / fdp;
    }

    return pstar;
}

double getustar(double pstar)
{
    double ustar;
    ustar = 0.5 * (ul + ur) + 0.5 * (f_r(pstar) - f_l(pstar));
    return ustar;
}

wavetype identifyNonlinearWave(double p)
{
    wavetype wave;
    if ( f(p) > 0)
    {
        wave = rarefaction;
    }
    else{
        wave = shock;
    }
    return wave;
}

double getrhostarl(double pstar, wavetype leftwave)
{
    double rhostarl;
    if (leftwave == rarefaction)
    {
        rhostarl = pstar / pl;
        rhostarl = rhol * pow(rhostarl,1.0/gamma);
    }
    else{
        rhostarl = pstar / pl + (gamma - 1) / (gamma + 1);
        rhostarl = rhostarl / ( ( (gamma-1)/(gamma+1) ) * ( pstar / pl) + 1 );
        rhostarl = rhol * rhostarl;
    }
    return rhostarl;
}

double getrhostarr(double pstar, wavetype rightwave)
{
    double rhostarr;
    if (rightwave == rarefaction)
    {
        rhostarr = pstar / pr;
        rhostarr = rhor * pow(rhostarr, 1.0/gamma);
    }
    else{
        rhostarr = pstar / pr + (gamma - 1) / (gamma + 1);
        rhostarr = rhostarr / ( ( (gamma-1)/(gamma+1) ) * (pstar / pr) + 1 );
        rhostarr = rhor * rhostarr;
    }
    return rhostarr;
}

void plot(wavetype leftwave, wavetype rightwave, double dll, double dlr, double dmid, double drl, double drr, double rhostarl, double rhostarr, double ustar, double pstar)
{
    int n = 1000;
    double h = (b-a) / n;
    double x[n], rhot, ut, pt, et;
    double locll, loclr, locmid, locrl, locrr;
    for (int i=0; i<n; i++)
    {
        x[i] = a +i*h;
    }

    locll = InitialDiscontinuity + dll;
    loclr = InitialDiscontinuity + dlr;
    locmid = InitialDiscontinuity + dmid;
    locrl = InitialDiscontinuity + drl;
    locrr = InitialDiscontinuity + drr;

    const char* frhon = "F:\\C++code\\RiemannSolver\\ExactRS\\ExactRSExactRS-density.plt";
    const char* fun = "F:\\C++code\\RiemannSolver\\ExactRS\\ExactRSExactRS-velocity.plt";
    const char* fpn = "F:\\C++code\\RiemannSolver\\ExactRS\\ExactRSExactRS-pressure.plt";
    const char* fen = "F:\\C++code\\RiemannSolver\\ExactRS\\ExactRSExactRS-internal_energy.plt";
    fstream fr;
    fstream fu;
    fstream fp;
    fstream fe;
    
    cout<<"locll="<<locll<<endl;
    cout<<"loclr="<<loclr<<endl;
    cout<<"locmid="<<locmid<<endl;
    cout<<"locrl="<<locrl<<endl;
    cout<<"locrr="<<locrr<<endl;
    remove(frhon);
    remove(fun);
    remove(fpn);
    remove(fen);
    fr.open(frhon, ios::out | ios :: app);
    fu.open(fun, ios::out | ios :: app);
    fp.open(fpn, ios::out | ios :: app);
    fe.open(fen, ios::out | ios :: app);
    fr<<"VARIABLES="<<"X"<<","<<"density"<<endl;
    fu<<"VARIABLES="<<"X"<<","<<"velocity"<<endl;
    fp<<"VARIABLES="<<"X"<<","<<"p"<<endl;
    fe<<"VARIABLES="<<"X"<<","<<"internal_energy"<<endl;

    for (int i=0; i<n; i++)
    {
        if (x[i] <= locll)
        {
            rhot = rhol;
            ut = ul;
            pt = pl;
            et = pt / ((gamma -1) * rhot);
        }
        else if (x[i] <= loclr && leftwave == rarefaction)
        {
            rhot = 2/(gamma + 1) + (gamma-1)*(ul - xi(x[i])/T)/( (gamma+1) * al );
            rhot = rhol * pow(rhot,2/(gamma-1));
            ut = 2 * (al + (gamma-1) * ul /2 + xi(x[i])/T) /(gamma + 1);
            pt = 2/(gamma + 1) + (gamma-1)*(ul - xi(x[i])/T) /( (gamma + 1) * al );
            pt = pl * pow(pt,2*gamma/(gamma-1));
            et = pt / ((gamma -1) * rhot);
        }
        else if (x[i] <= locmid)
        {
            rhot = rhostarl;
            ut = ustar;
            pt = pstar;
            et = pt / ((gamma -1) * rhot);
        }
        else if (x[i] <= locrl)
        {
            rhot = rhostarr;
            ut = ustar;
            pt = pstar;
            et = pt / ((gamma -1) * rhot);
        }
        else if (x[i] <= locrr && rightwave == rarefaction)
        {
            rhot = 2/(gamma+1) - (gamma-1)*(ur - xi(x[i])/T)/( (gamma+1)*ar );
            rhot = rhor * pow(rhot, 2/(gamma-1));
            ut = 2 * (-ar + (gamma - 1) * ur / 2 + xi(x[i])/T) / (gamma + 1);
            pt = 2/(gamma+1) - (gamma-1)*(ur - xi(x[i])/T)/( (gamma+1)*ar );
            pt = pr * pow(pt, 2*gamma/(gamma-1));
            et = pt / ((gamma -1) * rhot);
        }
        else
        {
            rhot = rhor;
            ut = ur;
            pt = pr;
            et = pt / ((gamma -1) * rhot);
        }
        fr<<"\t"<<x[i]<<"\t"<<rhot<<endl;
        fu<<"\t"<<x[i]<<"\t"<<ut<<endl;
        fp<<"\t"<<x[i]<<"\t"<<pt<<endl;
        fe<<"\t"<<x[i]<<"\t"<<et<<endl;
    }
    fr.close();
    fu.close();
    fp.close();
    fe.close();
    return ;
}

void exactRS()
{
    wavetype lwave, rwave;
    double pstar, ustar, rhostarl, rhostarr;
    double dll, dlr, dmid, drl, drr;
    
    pstar = getpstar();
    ustar = getustar(pstar);

    lwave = identifyNonlinearWave(pl);
    rwave = identifyNonlinearWave(pr);
    rhostarl = getrhostarl(pstar,lwave);
    rhostarr = getrhostarr(pstar,rwave);

    cout<<"leftwave="<<lwave<<endl;
    cout<<"rightwave="<<rwave<<endl;
    cout<<"0 for rarefaction wave, 1 for shock wave"<<endl<<endl;
    cout<<"pstar="<<pstar<<endl;
    cout<<"ustar="<<ustar<<endl;
    cout<<"rhostarl="<<rhostarl<<endl;
    cout<<"rhostarr="<<rhostarr<<endl<<endl;
    
    if (lwave == shock)
    {
        double s;
        s = (gamma + 1) * pstar / (2 * gamma * pl) + (gamma - 1)/(2 * gamma);
        s = ul - al * sqrt(s);

        dll = s * T;
        dlr = s * T;
    }
    else{
        double astarl, shl, stl;
        astarl = al * pow(pstar/pl,(gamma-1)/(2*gamma));

        shl = ul - al;
        stl = ustar - astarl;

        dll = shl * T;
        dlr = stl * T;
    }

    dmid = ustar * T;

    if (rwave == shock)
    {
        double s;
        s = (gamma + 1) * pstar/(2 * gamma * pr) + (gamma - 1)/(2 * gamma);
        s = ur + ar * sqrt(s);

        drl = s * T;
        drr = s * T;
    }
    else{
        double astarr, shr, str;
        astarr = ar * pow(pstar/pr, (gamma-1)/(2*gamma));
        shr = ur + ar;
        str = ustar + astarr;

        drl = str * T;
        drr = shr * T;
    }

    plot(lwave,rwave,dll,dlr,dmid,drl,drr,rhostarl,rhostarr,ustar,pstar);

    return ;
}

double* sampling(double x)
{
    double* ans = new double[4];
    wavetype lwave, rwave;
    double pstar, ustar, rhostarl, rhostarr;
    double rhot, ut, pt, et;
    double dll, dlr, dmid, drl, drr;
    double locll, loclr, locmid, locrl, locrr;
    
    pstar = getpstar();
    ustar = getustar(pstar);

    lwave = identifyNonlinearWave(pl);
    rwave = identifyNonlinearWave(pr);
    rhostarl = getrhostarl(pstar,lwave);
    rhostarr = getrhostarr(pstar,rwave);
    
    if (lwave == shock)
    {
        double s;
        s = (gamma + 1) * pstar / (2 * gamma * pl) + (gamma - 1)/(2 * gamma);
        s = ul - al * sqrt(s);

        dll = s * T;
        dlr = s * T;
    }
    else{
        double astarl, shl, stl;
        astarl = al * pow(pstar/pl,(gamma-1)/(2*gamma));

        shl = ul - al;
        stl = ustar - astarl;

        dll = shl * T;
        dlr = stl * T;
    }

    dmid = ustar * T;

    if (rwave == shock)
    {
        double s;
        s = (gamma + 1) * pstar/(2 * gamma * pr) + (gamma - 1)/(2 * gamma);
        s = ur + ar * sqrt(s);

        drl = s * T;
        drr = s * T;
    }
    else{
        double astarr, shr, str;
        astarr = ar * pow(pstar/pr, (gamma-1)/(2*gamma));
        shr = ur + ar;
        str = ustar + astarr;

        drl = str * T;
        drr = shr * T;
    }
    
    locll = InitialDiscontinuity + dll;
    loclr = InitialDiscontinuity + dlr;
    locmid = InitialDiscontinuity + dmid;
    locrl = InitialDiscontinuity + drl;
    locrr = InitialDiscontinuity + drr;


    if (x <= locll)
    {
        rhot = rhol;
        ut = ul;
        pt = pl;
        et = pt / ((gamma -1) * rhot);
    }
    else if (x <= loclr && lwave == rarefaction)
    {
        rhot = 2/(gamma + 1) + (gamma-1)*(ul - xi(x)/T)/( (gamma+1) * al );
        rhot = rhol * pow(rhot,2/(gamma-1));
        ut = 2 * (al + (gamma-1) * ul /2 + xi(x)/T) /(gamma + 1);
        pt = 2/(gamma + 1) + (gamma-1)*(ul - xi(x)/T) /( (gamma + 1) * al );
        pt = pl * pow(pt,2*gamma/(gamma-1));
        et = pt / ((gamma -1) * rhot);
    }
    else if (x <= locmid)
    {
        rhot = rhostarl;
        ut = ustar;
        pt = pstar;
        et = pt / ((gamma -1) * rhot);
    }
    else if (x <= locrl)
    {
        rhot = rhostarr;
        ut = ustar;
        pt = pstar;
        et = pt / ((gamma -1) * rhot);
    }
    else if (x <= locrr && rwave == rarefaction)
    {
        rhot = 2/(gamma+1) - (gamma-1)*(ur - xi(x)/T)/( (gamma+1)*ar );
        rhot = rhor * pow(rhot, 2/(gamma-1));
        ut = 2 * (-ar + (gamma - 1) * ur / 2 + xi(x)/T) / (gamma + 1);
        pt = 2/(gamma+1) - (gamma-1)*(ur - xi(x)/T)/( (gamma+1)*ar );
        pt = pr * pow(pt, 2*gamma/(gamma-1));
        et = pt / ((gamma -1) * rhot);
    }
    else
    {
        rhot = rhor;
        ut = ur;
        pt = pr;
        et = pt / ((gamma -1) * rhot);
    }

    ans[0] = rhot;
    ans[1] = ut;
    ans[2] = pt;
    ans[3] = et;

    return ans;
}