/**
 * @file sampling.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief Smpling states with exact Riemann solver
 * @version 0.1
 * @date 2022-05-21
 * 
 * @copyright Copyright (c) 2022 R.Y. Han, All Rights Reserved.
 * 
 */

#ifndef SAMPLING_H
#define SAMPLING_H

#include "configration.h"

/** @enum wavetype
 * @brief Define the notation of nonlinear waves.
 * 
 */
enum wavetype
{
    rarefaction = 0,    /**< rarefaction 0 */
    shock = 1           /**< shock 1 */
};

/**
 * @brief The function to compute the relation of pressure across the left nonlinear wave.
 * 
 * @param p : Pressure input.  
 * @param statesl : States of the left region in the local Riemann problem. Statesl = {rhol, ul, pl}
 * @return double : \n 
 * The input pressure p is larger than the pressure pl in the left region,
 * then use the function of a shock, else use the function of a rarefaction. Return the computed value.
 */
double f_l(double p, double* statesl);

/**
 * @brief The function to compute the relation of pressure across the right nonlinear wave.
 * 
 * @param p : Pressure input
 * @param statesr : States in the right region of the local Riemann problem. Statesr = {rhor, ur, pr}
 * @return double : \n 
 * If the input pressure p is larger than the pressure pl in the right region,
 * then use the function of a shock, else use the function of a rarefaction. Return the computed value.
 */
double f_r(double p, double* statesr);

/**
 * @brief f(p) = f_l(p) + f_r(p) + [u]
 * 
 * @param p : Pressure input.
 * @param statesl : Left states of the local Riemann problem
 * @param statesr : Right states of the local Riemann problem
 * @return double : \n 
 * When p = p*, f(p*) = 0.
 */
double f(double p, double* statesl, double* statesr);

/**
 * @brief The derivative function of f_l
 * 
 * @param p : Pressure input 
 * @param statesl : Left states of the local Riemann problem
 * @return double 
 */
double f_ldp(double p, double* statesl);

/**
 * @brief The derivative function of f_r
 * 
 * @param p : Pressure input 
 * @param statesr : Right states of the local Riemann problem
 * @return double 
 */
double f_rdp(double p, double* statesr);

/**
 * @brief Compute the pressure in the star region using Newton's iterative methods.
 * p* satisfies f(p*) = 0.
 * @param statesl : Left states of the local Riemann problem
 * @param statesr : Right states of the local Riemann problem
 * @return [double] p* : The pressure in the star region.
 */
double getpstar(double* statesl, double* statesr);

/**
 * @brief Compute the velocity in the star region.
 * 
 * @param pstar : the pressure in the star region 
 * @param statesl : Left states of the local Riemann problem
 * @param statesr : Right states of the local Riemann problem
 * @return u* : the velocity in the star region
 * @warning During iteration the parameter temp may be negative, choose the suitable initial value
 * of TRRS or mean value to avoid error.
 */
double getustar(double pstar, double* statesl, double* statesr);

/**
 * @brief Identify the type of two nonlinear waves.
 * @details 
 * -# If p* > pl, the left wave is a shock;
 * -# If p* < pl, the left wave is a rarefaction;
 * -# If p* > pr, the right wave is a shock;
 * -# If p* < pr, the right wave is a rarefaction.
 * 
 * @param p : pressure input, p = pl or pr
 * @param statesl : Left states of the local Riemann problem
 * @param statesr : Right states of the local Riemann problem
 * @return wavetype : 
 * -# When p = pl, return the type of left wave;
 * -# When p = pr, return the type of right wave.
 */
wavetype identifyNonlinearWave(double p, double* statesl, double* statesr);

/**
 * @brief Compute the density between left wave and the contact discontinuity.
 * 
 * @param pstar : the pressure in the star region
 * @param statesl : Left states of the local Riemann problem
 * @param statesr : Right states of the local Riemann problem
 * @param leftwave : the type of the left wave
 * @return double : the density between the left wave and the contact discontinuity.
 */
double getrhostarl(double pstar, double* statesl, wavetype leftwave);

/**
 * @brief Compute the density between right wave and the contact discontinuity.
 * 
 * @param pstar  : the pressure in the star region
 * @param statesl : Left states of the local Riemann problem
 * @param statesr : Right states of the local Riemann problem
 * @param rightwave : the type of the right wave
 * @return double : the density between the right wave and the contact discontinuity.
 */
double getrhostarr(double pstar, double* statesr, wavetype rightwave);

/**
 * @brief Get the required states at the given location x at the final time T.
 * 
 * @param x : given location
 * @param t : given time
 * @param statesl : Left states of the local Riemann problem
 * @param statesr : Right states of the local Riemann problem
 * @return The states at (x,T), a vector of 3 components: density, velocity, pressure.
 * @note 
 * -# In Godunov scheme we only need the solution at x=0.
 * -# Results are primary variables.
 * -# Given x = the considering location - the location of the discontinuity
 */
double* sampling(double x, double t, double* statesl, double* statesr);


#endif