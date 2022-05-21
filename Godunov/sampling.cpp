/**
 * @file sampling.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief Exact Riemann solver for the sampling process.
 * @version 1.01
 * @date 2022-05-21
 * 
 * @copyright Copyright (c) 2022 R.Y. Han, All Rights Reserved.
 * 
 */

#include "sampling.h"
#include <cmath>
#include <iostream>
using namespace std;

double f_l(double p, double* statesl)
{
    /**< statesl = {rhol, ul, pl} */
    double ans;
    double pl, rhol, Al, Bl;
    pl = statesl[2];
    rhol = statesl[0];
    Al = 2 / ( (gamma+1) * rhol);
    Bl = (gamma - 1) * pl / (gamma + 1);
    if (p > pl)
    {
        ans = Al / (p + Bl);
        ans = sqrt(ans);
        ans = ans * (p - pl);
    }
    else{
        double temp, al;
        temp = (gamma - 1)/(2 * gamma);
        rhol = statesl[0];
        al = sqrt( gamma * pl / rhol);

        ans = p / pl;
        ans = pow(ans, temp);
        ans = 2 * al * (ans - 1)/(gamma - 1);
    }
    return ans;
}

double f_r(double p, double* statesr)
{
    /** statesr = { rhor, ur, pr} */
    double ans;
    double pr, rhor, Ar, Br;
    pr = statesr[2];
    rhor = statesr[0];
    Ar = 2 / ((gamma + 1) * rhor);
    Br = (gamma - 1) * pr / (gamma + 1);
    if ( p > pr)
    {
        ans = Ar / ( p + Br );
        ans = sqrt(ans);
        ans = ans * (p - pr);
    }
    else{
        double temp, ar;
        temp = (gamma - 1) / (2 * gamma);
        ar = sqrt( gamma * pr / rhor);

        ans = p / pr;
        ans = pow(ans, temp);
        ans = 2 * ar * (ans - 1)/ (gamma - 1);
    }
    return ans;
}

double f(double p, double* statesl, double* statesr)
{
    double ans;
    ans = f_l(p,statesl) + f_r(p,statesr) + (statesr[1] - statesl[1]);  /** ur = statesl[1], ur = statesr[1] */
    return ans;
}

double f_ldp(double p, double* statesl)
{
    double ans;
    double pl, rhol, Al, Bl;
    pl = statesl[2];
    rhol = statesl[0];
    Al = 2 / ( (gamma+1) * rhol);
    Bl = (gamma - 1) * pl / (gamma + 1);
    if (p > pl)
    {
        ans = 1 - (p - pl)/(2 * (Bl + p));
        ans = ans * sqrt(Al/(Bl + p));
    }
    else{
        double temp, al;
        temp = (gamma + 1)/(2 * gamma);
        al = sqrt(gamma * pl / rhol);

        ans = p/pl;
        ans = pow(ans,-temp);
        ans = ans / (rhol * al);
    }
    return ans;
}

double f_rdp(double p, double* statesr)
{
    double ans;
    double pr, rhor, Ar, Br;
    pr = statesr[2];
    rhor = statesr[0];
    Ar = 2 / ((gamma + 1) * rhor);
    Br = (gamma - 1) * pr / (gamma + 1);
    if (p > pr)
    {
        ans = 1 - (p - pr)/(2 * (Br + p));
        ans = ans * sqrt(Ar/(Br + p));
    }
    else{
        double temp, ar;
        temp = (gamma + 1)/(2 * gamma);
        ar = sqrt(gamma * pr / rhor);

        ans = p/pr;
        ans = pow(ans,-temp);
        ans = ans / (rhor * ar);
    }
    return ans;
}

double getpstar(double* statesl, double* statesr)
{
    double pstar, temp;
    double fp, fdp, e, time;
    double pl, pr;
    pl = statesl[2];
    pr = statesr[2];
    e = 1e-6;
    time = 1;
    /** Initial pstar with TRRS */
    /*
    double ul, ur, rhol, rhor, al, ar;
    ul = statesl[1];
    ur = statesr[1];
    rhol = statesl[0];
    rhor = statesr[0];
    al = sqrt(gamma * pl / rhol);
    ar = sqrt(gamma * pr / rhor);
    pstar = al+ar - 0.5*(gamma-1)*(ur-ul);
    pstar = pstar / (pow(al/pl,(gamma-1)/(2*gamma)) + pow(ar/pr,(gamma-1)/(2*gamma)));
    pstar = pow(pstar,2*gamma/(gamma-1));
    //*/

    //** Initial pstar with mean value */
    pstar = 0.5 * (pl + pr);
    temp = 0;
    while(abs(pstar - temp)>e)
    {
        time++;
        if (time > 1e6)
        {
            break;
        }

        if (pstar <=0)
        {
            pstar = 1e-5;
        }
        temp = pstar;
        fp = f(temp,statesl, statesr);
        fdp = f_ldp(temp,statesl) + f_rdp(temp,statesr);
        pstar = temp - fp / fdp;
        
    }
    return pstar;
}

double getustar(double pstar, double* statesl, double* statesr)
{
    double ustar;
    double ul, ur;
    ul = statesl[1];
    ur = statesr[1];
    ustar = 0.5 * (ul + ur) + 0.5 * (f_r(pstar,statesr) - f_l(pstar,statesl));
    return ustar;
}

wavetype identifyNonlinearWave(double p, double* statesl, double* statesr)
{
    wavetype wave;
    if ( f(p,statesl, statesr) > 0)
    {
        wave = rarefaction;
    }
    else{
        wave = shock;
    }
    return wave;
}

double getrhostarl(double pstar, double* statesl, wavetype leftwave)
{
    double rhostarl;
    double rhol, pl;
    rhol = statesl[0];
    pl = statesl[2];
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

double getrhostarr(double pstar, double* statesr, wavetype rightwave)
{
    double rhostarr;
    double rhor, pr;
    rhor = statesr[0];
    pr = statesr[2];
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

double* sampling(double x, double t, double* statesl, double* statesr)
{
    double* ans = new double[4];
    wavetype lwave, rwave;
    double pstar, ustar, rhostarl, rhostarr;
    double rhot, ut, pt;
    double rhol, rhor, ul, ur, pl, pr, al, ar;
    double dll, dlr, dmid, drl, drr;
 
    rhol = statesl[0];
    rhor = statesr[0];
    ul = statesl[1];
    ur = statesr[1];
    pl = statesl[2];
    pr = statesr[2];
    al = sqrt(gamma * pl / rhol);
    ar = sqrt(gamma * pr / rhor);   
    
    if (rhol == rhor )
    {
        ans[0] = rhol;
        ans[1] = ul;
        ans[2] = pl;
        return ans;
    }
    pstar = getpstar(statesl,statesr);
    ustar = getustar(pstar,statesl, statesr);

    lwave = identifyNonlinearWave(pl,statesl,statesr);
    rwave = identifyNonlinearWave(pr,statesl,statesr);
    rhostarl = getrhostarl(pstar,statesl,lwave);
    rhostarr = getrhostarr(pstar,statesr,rwave);

    
    if (lwave == shock)
    {
        double s;
        s = (gamma + 1) * pstar / (2 * gamma * pl) + (gamma - 1)/(2 * gamma);
        s = ul - al * sqrt(s);

        dll = s;
        dlr = s;
    }
    else{
        double astarl, shl, stl;
        astarl = al * pow(pstar/pl,(gamma-1)/(2*gamma));

        shl = ul - al;
        stl = ustar - astarl;

        dll = shl;
        dlr = stl;
    }

    dmid = ustar;

    if (rwave == shock)
    {
        double s;
        s = (gamma + 1) * pstar/(2 * gamma * pr) + (gamma - 1)/(2 * gamma);
        s = ur + ar * sqrt(s);

        drl = s;
        drr = s;
    }
    else{
        double astarr, shr, str;
        astarr = ar * pow(pstar/pr, (gamma-1)/(2*gamma));
        shr = ur + ar;
        str = ustar + astarr;

        drl = str;
        drr = shr;
    }
    

    if (x/t <= dll)
    {
        rhot = rhol;
        ut = ul;
        pt = pl;
    }
    else if (x/t <= dlr && lwave == rarefaction)
    {
        rhot = 2/(gamma + 1) + (gamma-1)*(ul - x/t)/( (gamma+1) * al );
        rhot = rhol * pow(rhot,2/(gamma-1));
        ut = 2 * (al + (gamma-1) * ul /2 + x/t) /(gamma + 1);
        pt = 2/(gamma + 1) + (gamma-1)*(ul - x/t) /( (gamma + 1) * al );
        pt = pl * pow(pt,2*gamma/(gamma-1));
    }
    else if (x/t <= dmid)
    {
        rhot = rhostarl;
        ut = ustar;
        pt = pstar;
    }
    else if (x/t <= drl)
    {
        rhot = rhostarr;
        ut = ustar;
        pt = pstar;
    }
    else if (x/t <= drr && rwave == rarefaction)
    {
        rhot = 2/(gamma+1) - (gamma-1)*(ur - x/t)/( (gamma+1)*ar );
        rhot = rhor * pow(rhot, 2/(gamma-1));
        ut = 2 * (-ar + (gamma - 1) * ur / 2 + x/t) / (gamma + 1);
        pt = 2/(gamma+1) - (gamma-1)*(ur - x/t)/( (gamma+1)*ar );
        pt = pr * pow(pt, 2*gamma/(gamma-1));
    }
    else
    {
        rhot = rhor;
        ut = ur;
        pt = pr;
    }

    ans[0] = rhot;
    ans[1] = ut;
    ans[2] = pt;
    return ans;
}