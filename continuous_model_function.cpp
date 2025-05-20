#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "continuous_model_function.h"
#include "omp.h"



double normalize_angle(double angle) {
    angle = fmod(angle, 2 * M_PI);  // Réduit l'angle dans la plage [-2π, 2π]
    if (angle < 0) {
        angle += 2 * M_PI;  // Assure que l'angle est positif
    }
    return angle;
}

double todegree(double angle) {
    return angle * 180. / M_PI / 60;
}

void init_grid(t_spinner_grid* grid, int nx, int ny, double L){
    grid->nx = nx;
    grid->ny = ny;
    grid->Nspin = nx * ny;
    grid->L = L;

    grid->spin = (t_spinner*)malloc(sizeof(t_spinner) * grid->Nspin);
    if (grid->spin == NULL) {
		fprintf(stderr, "init_grid : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

    double Lx = L * cos(60. * M_PI / 180.);
    double Ly = L * cos(30. * M_PI / 180.);

    for(int i = 0; i < grid->nx; i++){
        for(int j = 0; j < grid->ny; j++){
            int index = j * grid->nx + i;
            grid->spin[index].x = 2. * Lx * i + Lx * (j % 2);
            grid->spin[index].y = - Ly * j;

            grid->spin[index].moments[0] = 1;
            grid->spin[index].moments[1] = 1;
            grid->spin[index].moments[2] = -1;

            grid->spin[index].omega = 0;
            grid->spin[index].theta_stat = 0;

            // droite, au dessus � droite, au dessus � gauche, gauche, en dessous � gauche, en dessous � droite
            
             if (i != 0) {
                grid->spin[index].neigbour[3] = j * grid->nx + i - 1; // voisin de gauche
            }
	        else { grid->spin[index].neigbour[3]  = - 1; }

            if (j != 0)     //peut avoir des voisins au dessus
            {
                if (j % 2 == 1)     //peut avoir un voisin au dessu � gauche
                {
                    grid->spin[index].neigbour[2] = (j - 1) * grid->nx + i;
                }
                else if (j % 2 == 0 && i > 0) {
                    grid->spin[index].neigbour[2] = (j - 1) * grid->nx + i - 1;
                }
	    	    else {
		    	    grid->spin[index].neigbour[2] = -1;
		        }

                if (j % 2 == 0)     //peut avoir un voisin au dessu � droite
                {
                    grid->spin[index].neigbour[1] = (j - 1) * grid->nx + i;
                }
                else if (j % 2 == 1 && i < grid->nx - 1) {
                    grid->spin[index].neigbour[1] = (j - 1) * grid->nx + i + 1;
		        }
		        else { grid->spin[index].neigbour[1] = -1; }
            }
	        else {
		        grid->spin[index].neigbour[2] = -1;
		        grid->spin[index].neigbour[1] = -1;
	        }

            if (i != grid->nx - 1) {
                grid->spin[index].neigbour[0] = j * grid->nx + i + 1;    //voisin de droite
	        }
	        else { grid->spin[index].neigbour[0] = -1; }

            if (j != grid->ny - 1)     //peut avoir des voisins en dessous
            {
                if (j % 2 == 1)     //peut avoir un voisin en dessous � gauche
                {
			        grid->spin[index].neigbour[4] = (j + 1) * grid->nx + i; 
                }
                else if (j % 2 == 0 && i > 0) {
                    grid->spin[index].neigbour[4] = (j + 1) * grid->nx + i - 1; 
		        }
		        else { grid->spin[index].neigbour[4] = -1; }

                if (j % 2 == 0)     //peut avoir un voisin en dessous � droite
                {
                    grid->spin[index].neigbour[5] = (j + 1) * grid->nx + i;
                }
                else if (j % 2 == 1 && i < grid->nx - 1) {
                    grid->spin[index].neigbour[5] = (j + 1) * grid->nx + i + 1;
		        }
		        else { grid->spin[index].neigbour[5] = -1; }
	        } else {// droite, au dessus � droite, au dessus � gauche, gauche, en dessous � gauche, en dessous � droite
		        grid->spin[index].neigbour[4] = -1;
		        grid->spin[index].neigbour[5] = -1;
	        }       
        }
    }
}

t_spinner_grid*  init_grid_array(double* allspin, int n_total, int nx, int ny, double L){
    t_spinner_grid* grid = (t_spinner_grid*)realloc(grid, sizeof(t_spinner_grid) * n_total);
    for(int i = 0; i < n_total; i++){
        init_grid(&grid[i], nx, ny, L);
        for(int j = 0; j < nx * ny; j++) {
            grid[i].spin[j].theta = allspin[i * nx * ny + j]; 
        }
    }
    return grid;
}

void init_rand(t_spinner_grid* grid, unsigned int* seed){

    for(int i = 0; i < grid->Nspin; i++) {
        grid->spin[i].theta = (double)rand_r(seed) / (double)RAND_MAX * 2. * M_PI;
        grid->spin[i].theta_stat = grid->spin[i].theta;
    }
    //for(int i = 0; i < grid->Nspin; i++) grid->spin[i].theta = 0;
}

void print_spinner(t_spinner_grid* grid) {

    for(int i = 0; i < grid->ny; i++) {
        for(int j = 0; j < grid->nx; j++) {
            printf("%f\t", grid->spin[i * grid->nx + j].theta);
        }
        printf("\n");
    }
}

double V_local(t_spinner_grid* grid, int index) {
    double Vtot = 0.0;  
    const double cst = MU0 / (4. * M_PI);  
    const double M_PI3 = 2. * M_PI /3.;

        for (int k = 0; k < 6; k++) {  

            int j = grid->spin[index].neigbour[k];

            if(j == -1) continue;

            for (int k = 0; k < NUMBER_ARM; k++) {
                for (int l = 0; l < NUMBER_ARM; l++) {

                    double m1_x = cos(grid->spin[index].theta + M_PI3 * k );
                    double m1_y = sin(grid->spin[index].theta + M_PI3 * k );

                    double m2_x = cos(grid->spin[j].theta + M_PI3 * l );
                    double m2_y = sin(grid->spin[j].theta + M_PI3 * l );

                    double dxkl = grid->spin[index].x + R * m1_x - grid->spin[j].x - R * m2_x ;
                    double dykl = grid->spin[index].y + R * m1_y - grid->spin[j].y - R * m2_y ;

                    double Lijkl = std::sqrt(dxkl * dxkl + dykl * dykl);

                    double Vinter = cst * (m1_x * m2_x + m1_y * m2_y - 3. * (m1_x * dxkl + m1_y * dykl) * (m2_x * dxkl + m2_y * dykl)
                                     /(Lijkl * Lijkl)) / (Lijkl * Lijkl * Lijkl);

                    Vtot += Vinter * grid->spin[index].moments[k] * grid->spin[j].moments[l];

                    //printf("%d %d %f %f %d %d %f %f %f a %f %f %f %f t %f %f\n", index, j, Vtot / HREF, Vinter /HREF, k, l, Lijkl, dxkl, dykl, m1_x, m1_y, m2_x, m2_y, grid->spin[j].theta, grid->spin[index].theta);
                }
            }
        }
    

    return Vtot / HREF ;
}

double V(t_spinner_grid* grid) {
    double Vtot = 0.0;  
    const double cst = MU0 / (4. * M_PI);  
    const double M_PI3 = 2. * M_PI /3.;

    for (int i = 0; i < grid->Nspin; i++) {
        for (int k = 0; k < 6; k++) {  

            int j = grid->spin[i].neigbour[k];

            if(j == -1) continue;

            for (int k = 0; k < NUMBER_ARM; k++) {
                for (int l = 0; l < NUMBER_ARM; l++) {

                    double m1_x = cos(grid->spin[i].theta + M_PI3 * k );
                    double m1_y = sin(grid->spin[i].theta + M_PI3 * k );

                    double m2_x = cos(grid->spin[j].theta + M_PI3 * l );
                    double m2_y = sin(grid->spin[j].theta + M_PI3 * l );

                    double dxkl = grid->spin[i].x + R * m1_x - grid->spin[j].x - R * m2_x ;
                    double dykl = grid->spin[i].y + R * m1_y - grid->spin[j].y - R * m2_y ;

                    double Lijkl = std::sqrt(dxkl * dxkl + dykl * dykl);

                    double Vinter = cst * (m1_x * m2_x + m1_y * m2_y - 3. * (m1_x * dxkl + m1_y * dykl) * (m2_x * dxkl + m2_y * dykl)
                                     /(Lijkl * Lijkl)) / (Lijkl * Lijkl * Lijkl);

                    Vtot += Vinter * grid->spin[i].moments[k] * grid->spin[j].moments[l];
                }
            }
        }
    }

    return Vtot / HREF / 2.;
}

void dV(t_spinner_grid* grid, double* grad, double delta_theta) {
    
    const double cst = MU0 / (4. * M_PI);

    for (int n = 0; n < grid->Nspin; n++) {
        
        double Vp = 0.0;  
        double Vm = 0.0;
      
        const double M_PI3 = 2. * M_PI /3.;

        for (int g = 0; g < 6; g++) {  

            int i = grid->spin[n].neigbour[g];

            if(i == -1) continue;  

            for (int k = 0; k < NUMBER_ARM; k++) {
                for (int l = 0; l < NUMBER_ARM; l++) {

                    double m1_x = cos(grid->spin[n].theta + delta_theta + M_PI3 * k );
                    double m1_y = sin(grid->spin[n].theta + delta_theta + M_PI3 * k );

                    double m2_x = cos(grid->spin[i].theta + M_PI3 * l );
                    double m2_y = sin(grid->spin[i].theta + M_PI3 * l );

                    double dxkl = grid->spin[n].x + R * m1_x - grid->spin[i].x - R * m2_x ;
                    double dykl = grid->spin[n].y + R * m1_y - grid->spin[i].y - R * m2_y ;

                    double Lijkl = std::sqrt(dxkl * dxkl + dykl * dykl);

                    Vp += cst * (m1_x * m2_x + m1_y * m2_y - 3. * (m1_x * dxkl + m1_y * dykl) * (m2_x * dxkl + m2_y * dykl)
                                    /(Lijkl * Lijkl))/ (Lijkl * Lijkl * Lijkl) * grid->spin[i].moments[k] * grid->spin[n].moments[l];
 
                }
            }

            for (int k = 0; k < NUMBER_ARM; k++) {
                for (int l = 0; l < NUMBER_ARM; l++) {

                    double m1_x = cos(grid->spin[n].theta - delta_theta  + M_PI3 * k );
                    double m1_y = sin(grid->spin[n].theta - delta_theta  + M_PI3 * k );

                    double m2_x = cos(grid->spin[i].theta + M_PI3 * l );
                    double m2_y = sin(grid->spin[i].theta + M_PI3 * l );

                    double dxkl = grid->spin[n].x + R * m1_x - grid->spin[i].x - R * m2_x ;
                    double dykl = grid->spin[n].y + R * m1_y - grid->spin[i].y - R * m2_y ;

                    double Lijkl = std::sqrt(dxkl * dxkl + dykl * dykl);

                    Vm += cst * (m1_x * m2_x + m1_y * m2_y - 3. * (m1_x * dxkl + m1_y * dykl) * (m2_x * dxkl + m2_y * dykl)
                                    /(Lijkl * Lijkl))/ (Lijkl * Lijkl * Lijkl) * grid->spin[i].moments[k] * grid->spin[n].moments[l];
                }
            }
        }
        grad[n] = (Vp - Vm) /2. /delta_theta;   
    }
}

void annealing(t_spinner_grid* grid, double T0, double Tf, double lambda, int niter, unsigned int* seed)
{

    grid->spin[0].theta = 4.77;
    grid->spin[1].theta = 1.61;

    int Niter = grid->Nspin * niter;
    
    for (double T = T0; T > Tf; T *= lambda) 
    {
        for (int j = 0; j < niter; j++)
        {
                int i = rand_r(seed) % grid->Nspin;

                double old_theta = grid->spin[i].theta;
                double old_V = V_local(grid, i);  // Old energy
          
                grid->spin[i].theta += 2. * ((rand_r(seed) / (double)RAND_MAX) - 0.5) * M_PI / 180.0 * 5.;  // Small random angle change
                //grid->spin[i].theta += (rand_r(seed) % 2) ? -0.5 : 0.5;  // Small random angle change

                double new_V = V_local(grid, i);  // New energy

                double delta_E = new_V - old_V;
                
                if (delta_E >= 0) 
                {
                    // If the new energy is higher, accept it with a certain probability
                    double acceptance_prob = exp(-delta_E / T);
                    if ( !(rand_r(seed) / (double)RAND_MAX < acceptance_prob) ) grid->spin[i].theta = old_theta;  // Revert to the old state

                }
            
        } //printf("%f %f %f %f %f\n",V(grid), V_local(grid, 0), V_local(grid, 1),V_local(grid, 2), V_local(grid, 3));
    }

    for(int j = 0; j < grid->Nspin; j++) grid->spin[j].theta = normalize_angle(grid->spin[j].theta );
    //printf("%f %f\n" , normalize_angle( grid->spin[0].theta), normalize_angle(grid->spin[1].theta));
}


void steepest_descent(t_spinner_grid* grid, double delta_theta, double alpha, int max_iter, double tolerance)
{
    double* grad = (double*)malloc(sizeof(double) * grid->Nspin);
    double* u = (double*)malloc(sizeof(double) * grid->Nspin);

    for(int i = 0; i < grid->Nspin; i++) u[i] = 1;
    double beta = 0.0;
    //alpha = 5;


    for(int i = 0; i < max_iter; i++)
    {   
        bool stop = true;

        //printf("%f %f\n" , normalize_angle( grid->spin[0].theta), normalize_angle(grid->spin[1].theta));
        
        dV(grid, grad, delta_theta);
        for(int j = 0; j < grid->Nspin; j++){
            
            u[j] = beta * u[j] - alpha * grad[j];
            grid->spin[j].theta += u[j];
            if (abs(alpha * grad[j]) > tolerance ) stop = false;

        }
        
        if (stop){
            //printf("break\n");
            break;
        }
    }

    for(int j = 0; j < grid->Nspin; j++) grid->spin[j].theta = normalize_angle(grid->spin[j].theta );

    free(grad);
    free(u);
}


bool isesqual(t_spinner_grid* grid1, t_spinner_grid* grid2, double tolerance) {
    for(int i = 0; i < grid1->Nspin; i++){
        if (abs(grid1->spin[i].theta - grid2->spin[i].theta) > tolerance) return false;
    }
    return true;
}

double* n_meta(int nx, int ny, double L, int n_iterations, double delta_theta, double alpha, int max_iter,
                                     double tolerance, unsigned int seed, int &n) {  

    t_spinner_grid* grid= (t_spinner_grid*)malloc(sizeof(t_spinner_grid) * n_iterations);
    bool* equal = (bool*)malloc(sizeof(bool) * n_iterations);

    for(int i = 0; i < n_iterations; i++){
        init_grid(&grid[i], nx, ny, L);
        init_rand(&grid[i], &seed);
        //steepest_descent(&grid[i], delta_theta, alpha, max_iter, tolerance);
        annealing(&grid[i], 1, 0.001, 0.5, 1000, &seed);
        equal[i] = false;
    }

    for(int i = 0; i < n_iterations; i++){
        for(int j = i + 1; j < n_iterations; j++){
            
            if(!equal[j] && isesqual(&grid[i], &grid[j], M_PI / 180. * 20.)){
                equal[j] = true;
            }
        }
    }

    n = 0;
    for(int i = 0; i < n_iterations; i++){
        if(!equal[i]) n++;
    }

    int N = nx * ny;
    double* spin = (double*)malloc(sizeof(t_spinner_grid) * n * N);
    int k = 0;
    for(int j = 0; j < n_iterations; j++){
        if(!equal[j]){
            for(int i = 0; i < N; i++) spin[k * N + i] = grid[j].spin[i].theta;
            k++;
        }
        free(grid[j].spin);
    }

    free(grid);
    free(equal);

    return spin;
}

double combine_std(double*  means, double* std_devs, int* sizes, int world_rank) {


    double mean = 0;
    int size = 0;

    for (int i = 0; i < world_rank; i++) {
        mean += means[i] * sizes[i];
        size += sizes[i];
    }
    mean /= size;

    double std = 0;
    for (int i = 0; i < world_rank; i++) std += (sizes[i] - 1) * 
                    (std_devs[i] * std_devs[i] + means[i] * means[i]);

    // Retourner l'écart type combiné (racine carrée de la variance combinée)
    return std::sqrt(std  / (size - 1) - mean * mean);
}



void n_meta_fct_inter(int nx, int ny, double L, int n_iterations, double delta_theta, double alpha, int max_iter, double tolerance) {

    // Déclarer le nombre de processus et leur rang dans MPI
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    unsigned int seed = time(NULL) + rank;

    // Répartir le travail entre les processus MPI
    int iterations_per_process = n_iterations / size;
    if (0 != n_iterations % size && rank == size - 1) iterations_per_process += n_iterations % size;

    // Variables pour stocker les résultats de chaque processus
    double energy = 0;
    double energy_std = 0;
    int* alln = (int*)calloc(size, sizeof(int));
    int n = 0;
    int n_total = 0;
    double* spin = n_meta(nx, ny, L, iterations_per_process, delta_theta, alpha, max_iter, tolerance, seed, n);

    MPI_Gather(&n, 1, MPI_INT, alln, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for(int i = 0 ; i < size; i++) n_total += alln[i];

    int* recvcounts = (int*)malloc(sizeof(int) * size);  // Le nombre d'éléments que chaque processus envoie
    int* displs = (int*)calloc(size, sizeof(int));       // Le décalage des éléments dans allspin

    for(int i = 0; i < size; i++) recvcounts[i] = alln[i] * nx * ny;  
    for(int i = 1; i < size; i++) displs[i] = displs[i - 1] + recvcounts[i - 1];  
    
    double* allspin = NULL;
    if (rank == 0)  allspin = (double*)malloc(sizeof(double) * n_total * nx * ny);

    MPI_Gatherv(spin, n * nx * ny, MPI_DOUBLE, allspin, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(recvcounts);
    free(displs);
    free(spin);

    if (rank == 0) {
        t_spinner_grid* grid = init_grid_array(allspin, n_total, nx, ny, L);
        free(allspin);

        bool* equal = (bool*)calloc(n_iterations, sizeof(bool));
        int n_diff = 0;
        for(int i = 0; i < n_total; i++){
            for(int j = i + 1; j < n_total; j++){
            
                if(!equal[j] && isesqual(&grid[i], &grid[j], M_PI / 180. * 20.)){
                    equal[j] = true;
                }
            }
        }

        for(int i = 0; i < n_total; i++){
            if(!equal[i]){
                energy += V(&grid[i]); //printf("%f\t%f\n", grid[i].spin[0].theta, grid[i].spin[1].theta);
                n_diff++; //printf("%f\n", V(&grid[i]));
            }
        }

        energy /= n_diff;

        for(int i = 0; i < n_total; i++){
            if(!equal[i]) energy_std += (V(&grid[i]) - energy) * (V(&grid[i]) - energy);
        }

        energy_std = std::sqrt(energy_std / (n_diff - 1));

        printf("%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\n", nx, ny, L, n_iterations, n_diff, 
               (double)n_diff / (double)n_iterations, energy / nx / ny, energy_std / nx / ny );
        
        for(int i = 0; i < n_total; i++) free(grid[i].spin);
        free(equal);
        free(grid);
    }

    // Synchroniser les processus MPI avant de finir
    MPI_Barrier(MPI_COMM_WORLD);

}



void discret(t_spinner_grid* grid, double mean[3]){
    double mcenter = 0;
    double medge = 0;
    double mcorner = 0;
    

    for( int i = 1; i < grid->nx - 1; i++){
        for( int j = 1; j < grid->ny - 1; j++){
            double est = normalize_angle(grid->spin[i].theta) * 180. / M_PI / 60.;
            double diff = est - (int)est;
            if (diff > 0.5) diff = 1 - diff;
            diff *= 60.;
            mcenter += diff;
        }
    }

    for( int i = 0; i < grid->nx; i++){
        for( int j = 0; j < grid->ny; j++){
            if(i == 0 || i == grid->nx - 1 || j == 0 || j == grid->ny - 1){
                double est = normalize_angle(grid->spin[i].theta) * 180. / M_PI / 60.;
                double diff = est - (int)est;
                if (diff > 0.5) diff = 1 - diff;
                diff *= 60.;
                if((i==0 && j==0) || (i==0 && j==grid->ny-1) || (i==grid->nx-1 && j==0) || (i==grid->nx-1 && j==grid->ny-1) ) mcorner += diff;
                else medge += diff;
            }
        }
    }
   
    if ((grid->nx - 2) * (grid->ny - 2) != 0) mcenter /= (grid->nx - 2) * (grid->ny - 2);
    mcorner /= 2 * (grid->nx != 1) + 2 * (grid->ny != 1); 
    if (2 * (grid->nx - 2) + 2 * (grid->ny - 2) != 0) medge /= 2 * (grid->nx - 2) + 2 * (grid->ny - 2);

    mean[0] = mcenter;
    mean[1] = medge;
    mean[2] = mcorner;

}

void discret_estim(int nx, int ny, int Nmean, int p){

    unsigned int seed = time(NULL);

    double* mean_center = (double*)malloc(sizeof(double) * Nmean);
    double* mean_edge = (double*)malloc(sizeof(double) * Nmean);
    double* mean_corner = (double*)malloc(sizeof(double) * Nmean);

    #pragma omp parallel num_threads(p)
    {
        spinner_grid spin;
        init_grid(&spin, nx, ny, 0.0319);

        unsigned int seedp = seed + omp_get_thread_num();

        #pragma omp for
        for(int i = 0; i < Nmean ; i++){
            init_rand(&spin, &seedp);
            steepest_descent(&spin, 5 * M_PI / 180., 0.001, 50000000, 1e-8 * M_PI /180);
            double mean[3];
            discret(&spin, mean);

            mean_center[i] = mean[0];
            mean_edge[i] = mean[1];
            mean_corner[i] = mean[2];
        }
        
        free(spin.spin);
    }

    double meanall[3] = {0.0, 0.0, 0.0};
    double sdall[3] = {0.0, 0.0, 0.0};

    for(int i = 0; i < Nmean; i++){
        meanall[0] += mean_center[i];
        meanall[1] += mean_edge[i];
        meanall[2] += mean_corner[i];
    }

    meanall[0] /= (double)Nmean;
    meanall[1] /= (double)Nmean;
    meanall[2] /= (double)Nmean;

    for(int i = 0; i < Nmean; i++){
        sdall[0] += (meanall[0] - mean_center[i]) * (meanall[0] - mean_center[i]);
        sdall[1] += (meanall[1] - mean_edge[i]) * (meanall[2] - mean_edge[i]);
        sdall[2] += (meanall[2] - mean_corner[i]) * (meanall[1] - mean_corner[i]);
    }

    sdall[0] = std::sqrt(sdall[0] / (double)Nmean);
    sdall[1] = std::sqrt(sdall[2] / (double)Nmean);
    sdall[2] = std::sqrt(sdall[1] / (double)Nmean);

    printf("%f\t%f\t%f\t%f\t%f\t%f\n", meanall[0], sdall[0], meanall[1], sdall[1], meanall[2], sdall[2]);

    free(mean_center);
    free(mean_corner);
    free(mean_edge);
}

