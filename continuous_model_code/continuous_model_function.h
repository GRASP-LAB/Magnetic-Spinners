#pragma once

#ifndef UNTITLED8_CONTINOUS_MODEL_FUNCTION_H
#define UNTITLED8_CONTINOUS_MODEL_FUNCTION_H

#endif //UNTITLED8_CONTINOUS_MODEL_FUNCTION_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <omp.h>
#include <ctime>

#define M_PI 3.14159265358979323846

#define MU0 0.0000012566370614 //[SI]
#define HREF 0.0497552 // d-d interaction for m = 1 and L = 0.0319 in the case where the interaction is maximal [SI] with mu_0/4Pi for mu = 1
#define R 0.008  // distance between the center of the dipole and the center of the spinner [m]


#define NUMBER_ARM 3

/**
 * @file continuous_model_function.h
 * @brief Functions and structures for simulating and analyzing spinner systems using continuous models.
 *
 * This header includes:
 * - Definitions for spinner and spinner grid structures
 * - Grid generation and initialization routines
 * - Energy calculations and gradient descent methods
 * - Simulated annealing for local minima exploration
 * - Stability analysis of spinner configurations
 * - Discretization error estimation between continuous and discrete models
 * - Parallel computing support with MPI and OpenMP
 */

/**
 * @brief structure for the spinner
 */
typedef struct spinner {
    double x;
    double y;
    double moments[3] = {1., 1., -1.};
    int neigbour[6];
    double theta; // radian
    double omega; // angular velocity

    double theta_stat; // friction

} t_spinner;

/**
 * @brief structure for the spinner grid
 */
typedef struct spinner_grid {
    int Nspin;
    double L;
    int nx;
    int ny;
    t_spinner* spin;
} t_spinner_grid;

/**
 * @brief normalise anglular coordinate 
 * 
 * @param  angle [IN] angle in R
 * 
 * @return angle in [0, 2 Pi] 
 */
double normalize_angle(double angle);

/**
 * @brief convert radian to degree
 */
double todegree(double angle);

/**
 * @brief generate the spinner grid with triangular unit cell 
 * 
 * @param  grid [OUT] spinner grid
 * 
 * @param  nx [IN] number of spinner on x-axis
 * 
 * @param  ny [IN] number of spinner on y-axis
 * 
 * @param  L [IN] lattice parametter
 */
void init_grid(t_spinner_grid* grid, int nx, int ny, double L);

/**
 * @brief generate the spinner grid with triangular unit cell
 * 
 * @param  allspin [IN] array of angle for each spinner
 * 
 * @param  n_total [IN] number of grid of spinner
 * 
 * @param  nx [IN] number of spinner on x-axis
 * 
 * @param  ny [IN] number of spinner on y-axis
 * 
 * @param  L [IN] lattice parametter
 * 
 * @return  return the array of spinner grid
 */
t_spinner_grid*  init_grid_array(double* allspin, int n_total, int nx, int ny, double L);

/**
 * @brief generate random configuration 
 * 
 * @param  grid [INOUT] spinner grid
 * 
 * @param  seed [IN] seed for the random number generation
 */
void init_rand(t_spinner_grid* grid, unsigned int* seed);

/**
 * @brief print the spinner grid
 */
void print_spinner(t_spinner_grid* grid);

/**
 * @brief calculate the torque for a spinner
 * 
 * @param  grid [IN] lattice of spinner
 * 
 * @param  index [IN] index of the spinner
 * 
 * @return  return the interaction energy for the spinner index
 */
double V_local(t_spinner_grid* grid, int index);

/**
 * @brief perfom the potential energy 
 * 
 * @param  grid [IN] lattice of spinner 
 * 
 * @return  return the potential enrgy total
 */
double V(t_spinner_grid* grid);

/**
 * @brief perfom the gradien of potential energy
 * 
 * @param grid [IN] lattice of spinenr 
 * 
 * @param  grad [OUT] gradient for each spinner 
 * 
 * @param delta_theta [IN] The discretisation for the cacluation of gradient
 */
void dV(t_spinner_grid* grid, double* grad, double delta_theta);


/**
 * @brief perfom the annealing of a spinner grid
 * 
 * @param  grid [INOUT] initial connfiguration and a the end this the final configuration
 * 
 * @param  T0 [IN] initial temperature
 * 
 * @param  Tf [IN] final temperature
 * 
 * @param  lambda [IN] cooling rate
 * 
 * @param  niter [IN] number of iteration for each temperature
 * 
 * @param  seed [IN] seed for the random number generation
 */
void annealing(t_spinner_grid* grid, double T0, double Tf, double lambda, int niter, unsigned int* seed);

/**
 * @brief steepest descent, find a configuration that is a local energy minimum 
 * 
 * @param  grid [INOUT] initial connfiguration and a the end this the final configuration
 * 
 * @param  delta_theta [IN] angle diffrance for the discret derivative of energy
 * 
 * @param  alpha [IN] steepest descent paramter
 * 
 * @param  max_inter [IN] the maxmimu number of iteration to find the local minima
 * 
 * @param tolerance [IN] if every spinner have a angle shift below this value, the funtion is stop
 */
void steepest_descent(t_spinner_grid* grid, double delta_theta, double alpha, int max_iter, double tolerance);

/**
 * @brief check if two spinner grid are equal
 * 
 * @param  grid1 [IN] first spinner grid
 * 
 * @param  grid2 [IN] second spinner grid
 * 
 * @param  tolerance [IN] the tolerance for the angle difference
 * 
 * @return  return true if the two grid are equal
 */
bool isesqual(t_spinner_grid* grid1, t_spinner_grid* grid2, double tolerance);

/**
 * @brief return proportion of stbale state renodmly generate
 * 
 * @param  nx [IN] x-size of lattice
 * 
 * @param  ny [IN] y-size of lattice
 * 
 * @param  L [IN] lattice parametter
 * 
 * @param  n_iterations [IN] number of simulation
 * 
 * @param  delta_theta [IN] angle diffrance for the discret derivative of energy
 * 
 * @param  alpha [IN] steepest descent paramter
 * 
 * @param  max_iter [IN] the maxmimu number of iteration to find the local minima
 * 
 * @param  tolerance [IN] if every spinner have a angle shift below this value, the funtion is stop
 * 
 * @param  seed [IN] seed for the random number generation
 * 
 * @param  n [OUT] number of different stable state
 * 
 * @return  return the number of different stable state
 */
double* n_meta(int nx, int ny, double L, int n_iterations, double delta_theta, double alpha, int max_iter,
                                     double tolerance, unsigned int seed, int &n);

/** 
 * @brief combine the standard deviation of different simulation
 * 
 * @param  means [IN] array of mean for each simulation
 * 
 * @param  std_devs [IN] array of standard deviation for each simulation
 * 
 * @param  sizes [IN] array of number of spinner for each simulation
 * 
 * @param  world_rank [IN]number of  the process
 * 
 * @return  return the combined standard deviation
 */
double combine_std(double*  means, double* std_devs, int* sizes, int world_rank);


/**
 * @brief search for metasbale state with MPI paralzization
 * 
 * @param  nx [IN] x-size of lattice
 * 
 * @param  ny [IN] y-size of lattice
 * 
 * @param  L [IN] lattice parametter
 * 
 * @param  n_iterations [IN] number of simulation
 * 
 * @param  delta_theta [IN] angle diffrance for the discret derivative of energy
 * 
 * @param  alpha [IN] steepest descent paramter
 * 
 * @param  max_iter [IN] the maxmimu number of iteration to find the local minima
 * 
 * @param  tolerance [IN] if every spinner have a angle shift below this value, the funtion is stop
 * 
 * @return  printf the result of the simulation
 */
void n_meta_fct_inter(int nx, int ny, double L, int n_iterations, double delta_theta, double alpha, int max_iter, double tolerance);

/**
 * @brief perfom the mean difference between continu model and discret model for center, edge and corner spinner 
 * 
 * @param  grid [OUT] spinner grid
 * 
 * @param  mean [OUT] mean[0] the diffrence for center spinner,  mean[1] for edge spinner, mean[2] for corner spinner 
 */
void discret(t_spinner_grid* grid, double mean[3]);

/**
 * @brief perfom the average mean difference between continu model and discret model for center, edge and corner spinner 
 * 
 * @param  nx [IN] x-size of lattice
 * 
 * @param  ny [IN] y-size of lattice
 * 
 * @param  Nmean [IN] number of simualtion for the average
 * 
 * @param  p [IN] number of threats for OpenMP
 */
void discret_estim(int nx, int ny, int Nmean, int p);


