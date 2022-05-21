/**
 * @file configuration.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief Initial the states of the Euler equations. Here we test the Riemann problems.
 * @version 1.01
 * @date 2022-05-21
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <cmath>

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

/**
 * @name Configuration
 *
 * @brief Configuration. Define the final time, location of the initial discontinuity,
 * left and right ends of the whole region, mapping to compute the slope,
 *  and ratio of the specific heats.
 * @{
 */
#define T 0.2 /**< final time */
#define InitialDiscontinuity 0.5    /**< the location of the initial discontinuity */
#define a 0.0   /**< left end */
#define b 1.0   /**< right end */

#define gamma 1.4   /** ratio of the specific heats. */
/** @} */

/**
* @name Initial states in the left region
* @brief Define the initial density, velocity, pressure
* and compute the sound speed in the left region.
* @{
*/
#define rhol_ini 1.0    /**< density in the left region */
#define ul_ini 0.0  /**< velocity in the left region */
#define pl_ini 1.0 /**< pressure in the left region */
#define al_ini sqrt( gamma * pl_ini / rhol_ini) /**< sound speed in the left region */
/** @} */ 

/**
 * @name Initial states in the right region
 * @brief Define the initial density, velocity, pressure
 *  and compute the sound speed of the right region.
 * @{
 */
#define rhor_ini 0.125    /**< density in the left region */
#define ur_ini 0.0  /**< velocity in the left region */
#define pr_ini 0.1    /**< pressure in the left region */
#define ar_ini sqrt(gamma * pr_ini / rhor_ini)  /**< sound speed in the left region */
/** @} */

/**
 * @brief Define the parameters of the uniform decomposition. Time step dt to be chosen.
 * 
 */
const int N=500;   /**< Total numbers of the cells */
const double h=(b-a)/N;   /**< Spacial step length */
const double CFL=0.5; /**< CFL number */

#endif