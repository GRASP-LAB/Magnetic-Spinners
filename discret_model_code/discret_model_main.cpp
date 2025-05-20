#include <iostream>
#include <omp.h>
#include "discret_model_function.h"
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <cstdint>


int main() {

    clock_t debut = clock();

    
    clock_t fin = clock();
    double tempsEcoule = (double)(fin - debut) / CLOCKS_PER_SEC;
    std::cout << "Temps ecoule : " << tempsEcoule << " secondes" << std::endl;

   
    return 0;
}

