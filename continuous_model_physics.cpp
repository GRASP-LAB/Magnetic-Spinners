#include <iostream>
#include "continuous_model_physics.h"
#include "continuous_model_function.h"
#include <fstream>



void B_dipole(double r, double rx, double ry, double mx, double my, double &bx, double &by) {

    // Calculer le produit scalaire (m · r̂)
    double m_dot_rhat = mx * rx / r + my * ry / r;

    double factor = MU * MU0 / (4 * M_PI * r * r * r);
    
    bx = factor * (3 * m_dot_rhat * rx / r - mx);
    by = factor * (3 * m_dot_rhat * ry / r - my);
}


double calculate_torque(t_spinner_grid *grid, int i, double bx, double by) {

    const double M_PI3 = 2. * M_PI /3.;

    double torque = 0.0;
    
    for (int n = 0; n < 6; n++) {  

        int j = grid->spin[i].neigbour[n];

        if(j == -1) continue;

        for (int k = 0; k < NUMBER_ARM; k++) {
            for (int l = 0; l < NUMBER_ARM; l++) {

                double m1_x = cos(grid->spin[i].theta + M_PI3 * k );
                double m1_y = sin(grid->spin[i].theta + M_PI3 * k );

                double m2_x = cos(grid->spin[j].theta + M_PI3 * l );
                double m2_y = sin(grid->spin[j].theta + M_PI3 * l );

                double dxkl = grid->spin[i].x + R * m1_x - grid->spin[j].x - R * m2_x ;
                double dykl = grid->spin[i].y + R * m1_y - grid->spin[j].y - R * m2_y ;

                double Lijkl = sqrt(dxkl * dxkl + dykl * dykl);

                double bxd, byd;
                B_dipole(Lijkl, dxkl, dykl, m2_x, m2_y, bxd, byd);

                torque += MU * (m1_x * byd - m1_y * bxd) * grid->spin[j].moments[l] *grid->spin[i].moments[k];

            }
        }
    }

    for (int k = 0; k < NUMBER_ARM; k++) {
        torque += MU * (cos(grid->spin[i].theta + M_PI3 * k ) * by - sin(grid->spin[i].theta + M_PI3 * k ) * bx) * grid->spin[i].moments[k];
    }

    return torque;
}


void print_physics(t_spinner_grid *grid, t_params *params, double bx, double by) {

    FILE* fichier = fopen(params->filename, "w");

    if (fichier == nullptr) std::cout << "Impossible d'ouvrir le fichier." << std::endl;

    fprintf(fichier, "%d\t%d\t%f\n", grid->nx, grid->ny, grid->L);

    double t = 0;
    int index = 0;
    while(t < params->tmax) {

        for (int i = 0; i < grid->Nspin; i++) {
        
            double torque = calculate_torque(grid, i, bx, by);

            if(fabs(torque) < 0.000001) torque = 0;
            
            double alpha = (torque -  FRICTION * grid->spin[i].omega) / INERTIE; 
            grid->spin[i].omega += alpha * params->dt;
            
        }

        for (int i = 0; i < grid->Nspin; i++) {
            grid->spin[i].theta += grid->spin[i].omega * params->dt / 2.; 
        }

        t += params->dt;
        index++;

        if (index % params->rate_sample == 0 || index == 1) {    
            fprintf(fichier, "%f\t", t);
            for (int i = 0; i < grid->Nspin - 1; i++) {
                fprintf(fichier, "%f\t", normalize_angle(grid->spin[i].theta));
            }
            fprintf(fichier, "%f\n", normalize_angle(grid->spin[grid->Nspin - 1].theta));
        }
    }
    
    fclose(fichier);
}


