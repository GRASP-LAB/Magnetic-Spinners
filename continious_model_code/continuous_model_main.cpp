#include <iostream>
#include <mpi.h>
#include <omp.h>
#include "continuous_model_function.h"
#include "continuous_model_physics.h"
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <cstring>
#include "continuous_model_um.h"


int main() {

    clock_t debut = clock();

    MPI_Init(NULL, NULL);


    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

    MPI_Finalize();
   
    return 0;
}

