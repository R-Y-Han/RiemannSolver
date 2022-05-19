/**
 * @file configration.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief Configrations of the 1st order upwind schemes.
 * @version 1.01
 * @date 2022-05-17
 * 
 * @copyright Copyright (c) 2022 R.Y. Han, All Rights Reserved.
 * 
 */

#ifndef CONFIGRATION_H
#define CONFIGRATION_H

#define a 1.0     /**< characteristic speed */
#define left -1.0    /**< the left end of the initial interval */
#define right 1.0     /**< the right end of the initial interval */
#define n 200   /**< number of decomposition cells  */
#define c 0.8 /**< CFL number */
#define nt  125 /**< number of time steps */

const double h = (right - left)/n;  /**< spatial step   */
const double dt = h * c; /**< time step  */
const double T = dt * nt;    /**< finial time */

/**
 * @brief The initial value function.
 * 
 * @param x 
 * @return u_0(x)
 */
double u_0(double x);

/**
 * @brief The flux function
 * 
 * @param u 
 * @return f(u)
 */
double f(double u);

/**
 * @brief conpute the exact solution to the equation u_t + u_x = 0
 * 
 * @param x 
 * @param t 
 * @return u(x,t)
 * @note for linear equation
 */
double u_exact(double x, double t);

#endif