#pragma once

#ifndef UNTITLED8_CONTINOUS_MODEL_PHYSICS_H
#define UNTITLED8_CONTINOUS_MODEL_PHYSICS_H

#endif //UNTITLED8_CONTINOUS_MODEL_PHYSICS_H

#include "continuous_model_function.h"

#define MU 0.02278125 // moment magnetique d'un aimant
#define INERTIE 1e-8 // moment of inertia
#define FRICTION 0.00005 * 0.05 // roulement a bille

/**
 * @brief structure for the parameters
 */
typedef struct params
{
    double dt;
    double tmax;
    char *filename;
    int rate_sample;

}t_params;

void B_dipole(double r, double rx, double ry, double mx, double my, double &bx, double &by);

double calculate_torque(t_spinner_grid *grid, int i, double bx, double by);

void print_physics(t_spinner_grid *grid, t_params *params, double bx, double by);