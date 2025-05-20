#pragma once

#ifndef UNTITLED8_CONTINOUS_MODEL_UM_H
#define UNTITLED8_CONTINOUS_MODEL_UM_H

#endif // UNTITLED8_CONTINOUS_MODEL_UM_H

#include <iostream>
#include <climits>
#include <float.h>
#include "continuous_model_function.h"
#include <cstring>




typedef struct matrice_line{
	double* col;
	int pos;
}matrice_line_t;

typedef struct matrice{
	int N;
	matrice_line_t* line;
}matrice_t;

typedef struct cluster{
	int size;
	int* pos;
	bool notmerged;
}cluster_t;

typedef struct tree{
	int N;
	cluster_t* cluster;
}tree_t;

t_spinner_grid* generate_state( int nx, int ny, double L, unsigned int* seed);

double dist_EL(t_spinner_grid* grid, int i, int j);

double single_link(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist );

void cluster_fusion(tree_t* tree, matrice_t* matrice_dist, matrice_t* matrice_ultra, double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ));

void matrice_ultra(matrice_t* matrice_dist, matrice_t* matrice_ultra, double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ));

double degree_UM(t_spinner_grid* grid, int Ngrid, double (*dist)(t_spinner_grid* grid, int i, int j), 
				double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ) );

void print_matrice(matrice_t* matrice, FILE* fichier);

FILE* openfile_out(char* add) ;

double distline(double* A, double* B, int N);

void tri(matrice_t* matrice, matrice_t* ultra);

void print_dist(t_spinner_grid* grid, int Ngrid, char* add, double (*dist)(t_spinner_grid* grid, int i, int j), 
		double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ));