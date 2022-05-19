/**
 * @file configration.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief The detail definitions of configration.h
 * @version 1.01
 * @date 2022-05-17
 * 
 * @copyright Copyright (c) 2022 R.Y. Han, All Rights Reserved.
 * 
 */

#include "configration.h"
#include <cmath>

double u_0(double x)
{
    double ans;
    ans = exp(-x*x);
    return ans;
}

double f(double u)
{
    return a*u;
}

double u_exact(double x, double t)
{
    return u_0(x-a*t);
}