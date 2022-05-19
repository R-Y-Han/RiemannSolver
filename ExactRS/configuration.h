/**
 * @file configuration.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief Configuration of the Riemann problem.
 * @version 1.01
 * @date 2022-05-14
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
#define T 0.035 /**< final time */
#define InitialDiscontinuity 0.5    /**< the location of the initial discontinuity */
#define a 0.0   /**< left end */
#define b 1.0   /**< right end */
#define xi(x) ( (x - InitialDiscontinuity) )    /**< mapping */

#define gamma 1.4   /** ratio of the specific heats. */
/** @} */

/**
* @name States in the left region
* @brief Define the initial density, velocity, pressure
* and compute the sound speed in the left region.
* @{
*/
#define rhol 1.0    /**< density in the left region */
#define ul 0.0  /**< velocity in the left region */
#define pl 0.01 /**< pressure in the left region */
#define al sqrt( gamma * pl / rhol) /**< sound speed in the left region */
/** @} */ 

/**
 * @name States in the right region
 * @brief Define the initial density, velocity, pressure
 *  and compute the sound speed of the right region.
 * @{
 */
#define rhor 1.0    /**< density in the left region */
#define ur 0.0  /**< velocity in the left region */
#define pr 100.0    /**< pressure in the left region */
#define ar sqrt(gamma * pr / rhor)  /**< sound speed in the left region */
/** @} */

/**
 * @brief Define the coefficients used in the computation of function f_l, f_r,
 * @see :: f_l, f_r
 * @name Some coefficients
 */
#define Al 2 / ( (gamma+1) * rhol)
#define Bl (gamma - 1) * pl / (gamma + 1)
#define Ar  2 / ((gamma + 1) * rhor)
#define Br  (gamma - 1) * pr / (gamma + 1)
/** @} */

/** @enum wavetype
 * @brief Define the notation of nonlinear waves.
 * 
 */
enum wavetype
{
    rarefaction = 0,    /**< rarefaction 0 */
    shock = 1           /**< shock 1 */
};


#endif