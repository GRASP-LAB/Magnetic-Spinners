#pragma once

#ifndef UNTITLED8_CONTINOUS_MODEL_UM_H
#define UNTITLED8_CONTINOUS_MODEL_UM_H

#endif // UNTITLED8_CONTINOUS_MODEL_UM_H

#include <iostream>
#include <climits>
#include <float.h>
#include "continuous_model_function.h"
#include <cstring>


/**
 * @file cluster_analysis.h
 * @brief Functions related to spinner grid generation, clustering, and ultrametric distance computation.
 *
 * This header declares functions for:
 * - Generating a grid of magnetic spinners
 * - Computing distances and similarities between elements
 * - Performing hierarchical clustering (e.g., single-link)
 * - Building an ultrametric matrix from a distance matrix
 * - Measuring the ultrametric degree
 * - Utilities for file I/O and matrix manipulation
 */

/**
 * @brief structure for the distance matrix line
 */
typedef struct matrice_line{
	double* col;
	int pos;
}matrice_line_t;

/**
 * @brief structure for the distance matrix
 */
typedef struct matrice{
	int N;
	matrice_line_t* line;
}matrice_t;

/**
 * @brief structure for the cluster for clustering
 */
typedef struct cluster{
	int size;
	int* pos;
	bool notmerged;
}cluster_t;

/**
 * @brief structure for the tree of cluster
 */
typedef struct tree{
	int N;
	cluster_t* cluster;
}tree_t;


/**
 * @brief Generate a spinner grid state.
 * 
 * @param nx Number of spinners in x direction.
 * @param ny Number of spinners in y direction.
 * @param L Lattice spacing.
 * @param seed Pointer to a random seed.
 * @return Pointer to the generated spinner grid.
 */
t_spinner_grid* generate_state(int nx, int ny, double L, unsigned int* seed);

/**
 * @brief Compute the  distance based on difference of local energy between two elements in the spinner grid.
 * 
 * @param grid Pointer to the spinner grid.
 * @param i Index of first element.
 * @param j Index of second element.
 * @return Distance between element i and j.
 */
double dist_EL(t_spinner_grid* grid, int i, int j);

/**
 * @brief Compute the single-link similarity between two clusters.
 * 
 * @param cluster_A First cluster.
 * @param cluster_B Second cluster.
 * @param matrice_dist Distance matrix.
 * @return Similarity between the two clusters.
 */
double single_link(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist);

/**
 * @brief Merge two clusters in the tree using the given similarity function.
 * 
 * @param tree Pointer to the tree structure.
 * @param matrice_dist Distance matrix.
 * @param matrice_ultra Ultrametric matrix (output).
 * @param similarity Similarity function to use (e.g., single-link).
 */
void cluster_fusion(tree_t* tree, matrice_t* matrice_dist, matrice_t* matrice_ultra, double (similarity)(cluster_t, cluster_t, matrice_t*));

/**
 * @brief Compute the ultrametric matrix from a distance matrix.
 * 
 * @param matrice_dist Original distance matrix.
 * @param matrice_ultra Ultrametric matrix to be filled.
 * @param similarity Similarity function to use (e.g., single-link).
 */
void matrice_ultra(matrice_t* matrice_dist, matrice_t* matrice_ultra, double (similarity)(cluster_t, cluster_t, matrice_t*));

/**
 * @brief Calculate the ultrametric degree of the system.
 * 
 * @param grid Spinner grid.
 * @param Ngrid Number of elements in the grid.
 * @param dist Function to compute distance between elements.
 * @param similarity Similarity function between clusters.
 * @return Degree of ultrametricity.
 */
double degree_UM(t_spinner_grid* grid, int Ngrid, double (*dist)(t_spinner_grid*, int, int), 
                 double (similarity)(cluster_t, cluster_t, matrice_t*));

/**
 * @brief Print the contents of a matrix to a file.
 * 
 * @param matrice Matrix to print.
 * @param fichier Output file pointer.
 */
void print_matrice(matrice_t* matrice, FILE* fichier);

/**
 * @brief Open a file for output writing.
 * 
 * @param add Path to the file.
 * @return File pointer.
 */
FILE* openfile_out(char* add);

/**
 * @brief Compute the Euclidean distance between two vectors.
 * 
 * @param A First vector.
 * @param B Second vector.
 * @param N Size of vectors.
 * @return Distance between A and B.
 */
double distline(double* A, double* B, int N);

/**
 * @brief Sort or process two matrices, potentially aligning ultrametric and original matrices.
 * 
 * @param matrice Input matrix.
 * @param ultra Ultrametric matrix.
 */
void tri(matrice_t* matrice, matrice_t* ultra);

/**
 * @brief Print distance matrix and associated ultrametric results to a file.
 * 
 * @param grid Spinner grid.
 * @param Ngrid Number of elements in grid.
 * @param add Output file path.
 * @param dist Distance function.
 * @param similarity Similarity function.
 */
void print_dist(t_spinner_grid* grid, int Ngrid, char* add, 
                double (*dist)(t_spinner_grid*, int, int), 
                double (similarity)(cluster_t, cluster_t, matrice_t*));
