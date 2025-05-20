#pragma once

#ifndef UNTITLED8_DISCRET_MODEL_FUNCTION_H
#define UNTITLED8_DISCRET_MODEL_FUNCTION_H

#endif //UNTITLED8_DISCRET_MODEL_FUNCTION_H

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#define PI 3.14159265358979323846 
#define PI3 3.14159265358979323846/3. 
#define MU0 0.0000012566370614 //[SI]
#define HREF 0.0497552 // d-d interaction for m = 1 and L = 0.0319 in the case of maximal interaction [SI], using mu_0 / (4π) for mu = 1
#define H25MAX 1.258910 // maximum element of the H tensor for L = 25 mm

#define MU 0.02278125 // magnetic moment of a magnet

#define META_ANGLE 3.2 // angle in degrees used to perturb a spinner and test its stability

#define SIZE_H 6*6
#define SIZE_NEIGHBOUR 6
#define R 0.008 // distance between the center of the dipole and the center of the spinner [m]
#define STRING_MAX 256

#define DBL_PRECISION 1.0E-14 // precision threshold for metastability criterion

#define GETELEMENT(arr, i, v) arr[ i * 6 + v]
#define SETELEMENT(arr, i, v, value) arr[ i * 6 + v] = value

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

typedef struct spinners {
	
	int nx;
	int ny;
	int Ngrid;
	double L;
	int* angles;
	int* anti;

} spinners_t;

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





/********************************************************************************/
/*                                                                              */
/*                            comput Energy dipole fct                          */
/*                                                                              */
/********************************************************************************/

/**
 * @brief compute a dipole dipole interaction 
 * 
 * @param  m1_x [IN] composante x of dipole 1
 * 
 * @param  m1_y [IN] composante y of dipole 1
 * 
 * @param  m2_x [IN] composante x of dipole 2
 * 
 * @param  m2_y [IN] composante y of dipole 2
 * 
 * @param  r1_x [IN] composante x of position of dipole 1
 * 
 * @param  r1_y [IN] composante y of position of dipole 1
 * 
 * @param  r2_x [IN] composante x of position of dipole 2
 * 
 * @param  r2_y [IN] composante y of position of dipole 2
 * 
 * @return energy of dipole dipole interaction 
 */
double Udipole(double m1_x, double m1_y, double m2_x, double m2_y, double r1_x, double r1_y, double r2_x, double r2_y);


/**
 * @brief compute a dipole B-fiedl interaction
 *
 * @param  mx [IN] composante x of dipole
 *
 * @param  my [IN] composante y of dipole 
 *
 * @param  bx [IN] composante x of B-field
 *
 * @param  by [IN] composante y of B-field
 *
 * @return energy of dipole B-field interaction -mu.B
 */
double UB(double mx, double my, double bx, double by);

/**
 * @brief compute a ppm spinner spinner interaction
 *
 * @param  theta [IN] angle of spinner of interest
 *
 * @param  thetav [IN]  angle of niegtboor of spinner of interest
 *
 * @param  L [IN] distance between 2 spinners
 *
 * @param  nv [IN] position of neigtboor : right 1 to Below right 5 by counterclockwise
 *
 * @return energy of ppm spinner spinner interaction normalised by HERF
 */
double Uspinner(int theta, int thetav, double L, int nv);


/**
 * @brief compute a ppm spinner spinner interaction, for a +++ spinner
 *
 * @param  theta [IN] angle of spinner of interest
 *
 * @param  thetav [IN]  angle of niegtboor of spinner of interest
 *
 * @param  L [IN] distance between 2 spinners
 *
 * @param  nv [IN] position of neigtboor : right 1 to Below right 5 by counterclockwise
 *
 * @return energy of ppm spinner spinner interaction normalised by HERF
 */
double Uspinner_ppp(int theta, int thetav, double L, int nv);

double Uspinner_pertubed(int theta, int thetav, double L, int nv, double beta);


/**
 * @brief compute a ppm spinner B-field interaction
 *
 * @param  theta [IN] angle of spinner of interest
 *
 * @param  bx [IN] composante x of B-field
 *
 * @param  by [IN] composante y of B-field
 *
 * @return energy of ppm spinner B-field interaction normalised by HERF
 */
double UBspinner(int theta, double bx, double by);

/**
 * @brief compute the interaction tensor by dipole dipole interaction
 *
 * @param  L [IN] distance between 2 spinner
 *
 * @return interaction tensor of the right neghitbourg
 */
double* H_init(double L);

/**
 * @brief compute the interaction tensor by dipole dipole interaction. The erngy is renormisze by the -- interaction of the L considerded
 *
 * @param  L [IN] distance between 2 spinner
 *
 * @return interaction tensor of the right neghitbourg
 */
double*  H_init_L(double L);

double*  H_init_L25(double L);


/**
 * @brief compute the interaction tensor by dipole dipole interaction for a +++ spinener
 *
 * @param  L [IN] distance between 2 spinner
 *
 * @return interaction tensor of the right neghitbourg for +++ spinner
 */
double*  H_ppp_init(double L);

/**
 * @brief compute the interaction tensor  by dipole B-field interaction
 *
 * @param  bx [IN] composante x of B-field
 *
 * @param  by [IN] composante y of B-field
 *
 * @return interaction tensor of the right neghitbourg
 */
double* H_B_init(double bx, double by);

/**
 * @brief compute the interaction tensor  by dipole B-field interaction, it renormalized by HREF for the specific L
 *
 * @param  bx [IN] composante x of B-field
 *
 * @param  by [IN] composante y of B-field
 * 
 * @param  L [IN] lenght betwee two spinner center
 * 
 * @return interaction tensor of the right neghitbourg
 */
double* H_B_init_L( double bx, double by, double L);

/**
 * @brief plot interaction tensor 
 *
 * @param  H [IN]  interaction tensor dipole dipole
 */
void H_plot(double* H);


/**
 * @brief plot interaction tensor
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void H_B_plot(double* HB);


/********************************************************************************/
/*                                                                              */
/*          Initialisation spinner lattice and lattice fct                      */
/*                                                                              */
/********************************************************************************/

/**
 * @brief intitialize a spinner lattice structure
 *
 * @param  spin [OUT]  spinner_t
 * 
 * @param  L [IN]  distance between tow spinners
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 * 
 * @param  Ngrid [IN]  number de grid
 */
void spinners_init(spinners_t* spin, double L, int nx, int ny, int Ngrid);

/**
 * @brief intitialize a spinner lattice structure with spinner PPM and MPP
 *
 * @param  spin [OUT]  spinner_t
 * 
 * @param  L [IN]  distance between tow spinners
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 * 
 * @param  seed [IN]  seed for rand_d threat safe
 */
void spinners_init_anti(spinners_t *spin, double L, int nx, int ny, unsigned int* seed, int anti);

/**
 * @brief intitialize a spinner lattice structure with spinner PPM and MPP
 *
 * @param  spin [OUT]  spinner_t
 * 
 * @param  L [IN]  distance between tow spinners
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 * 
 * @param  seed [IN]  seed for rand_d threat safe
 * 
 * @param  prob_anti [IN]  probability that a spinner has ta have -1 as dipole
 */
void spinners_init_anti_prob(spinners_t *spin, double L, int nx, int ny, int Ngri, unsigned int* seed, double prob_anti);

/**
 * @brief free interaction tensor and spinners_t
 *
 * @param  spin [INOUT]  spinner_t
 * 
 * @param  H [INOUT]  interaction tensor
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void Finalisation_simu(spinners_t* spin, double* H, double* HB);

/**
 * @brief comput the index of spinner neighbour
 *
 * @param  spin [IN]  spinner lattice
 * 
 * @param  index [IN]  index of interesst spinner
 * 
 * @return array with neighbour position
 */
int* neighbour( spinners_t* spin, int index);


/********************************************************************************/
/*                                                                              */
/*                           magnetisation                                      */
/*                                                                              */
/********************************************************************************/

/**
 * @brief comput the magnetisation of a spinner lattice
 *
 * @param  spin [IN]  spinner lattice
 * 
 * @return magnetisation of the spinner lattice
 */
double magnetisation(spinners_t* spin);

/********************************************************************************/
/*                                                                              */
/*                           Energy of spinner lattice                          */
/*                                                                              */
/********************************************************************************/


/**
 * @brief comput the Energy betwenn a spinner and all neighbour
 *
 * @param  spin [IN]  spinner lattice
 *
 * @param  index [IN]  index of interesst spinner
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return local Energy
 */
double E_local(spinners_t* spin, int index, double* H, double* HB, int offset);

/**
 * @brief comput the Energy betwenn a spinner and all neighbour
 *
 * @param  spin [IN]  spinner lattice
 *
 * @param  index [IN]  index of interesst spinner
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return local Energy
 */
double E_local_anti(spinners_t* spin, int index, double* H, double* HB, int offset);


/**
 * @brief comput the total Energy of spinners grid
 *
 * @param  spin [IN]  spinner lattice
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return total Energy
 */
double E_total(spinners_t* spin, double* H, double* HB, int offset);

/**
 * @brief comput the total Energy of spinners grid
 *
 * @param  spin [IN]  spinner lattice
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return total Energy
 */
double E_total_anti(spinners_t* spin, double*H, double* HB, int offset);

/**
 * @brief detremine if spinner configuration is metastable
 * 
 * @param  spin [IN]  spinner lattice
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return true if spinner configuration is metastable, or false
 */
bool metastable(spinners_t* spin, double* H, double* HB, int offset);

/**
 * @brief detremine if spinner configuration is metastable
 * 
 * @param  spin [IN]  spinner lattice
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return true if spinner configuration is metastable, or false
 */
bool metastable_anti(spinners_t* spin, double* H, double* HB, int offset);

/********************************************************************************/
/*                                                                              */
/*         						 to continu     			                    */
/*                                                                              */
/********************************************************************************/

/**
 * @brief check if a spinner configuration is metastable
 *
 * @param  spin [IN]  spinner lattice
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return true if spinner configuration is metastable, or false
 */
bool meta_continu(spinners_t* spin, double* H, double* HB,  int offset);

/**
 * @brief check if a spinner configuration is metastable
 *
 * @param  spin [IN]  spinner lattice
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 *
 * @return true if spinner configuration is metastable, or false
 */
void to_meta_continu(spinners_t* spin);


/********************************************************************************/
/*                                                                              */
/*                                      utilitaire                              */
/*                                                                              */
/********************************************************************************/


/**
 * @brief interverted tow int
 *
 * @param  A [INOUT]   int that will be interverted 
 * 
 * @param  B [INOUT]   int that will be interverted 
 */
void change(int* A, int* B);

/**
 * @brief compute the mean of an array
 *
 * @param  A [IN]  array
 * 
 * @param  N [IN]  size of the array
 *
 * @return mean of the array
 */
double mean(double* A, int N);

/**
 * @brief compute the standard deviation of an array
 *
 * @param  A [IN]  array
 * 
 * @param  N [IN]  size of the array
 *
 * @return standard deviation of the array
 */
double SD(double* A, int N);

/**
 * @brief open a file in w write and erase mode
 *
 * @param  add [IN]  pathfile
 *
 * @return a FILE*
 */
FILE* openfile_out(char* add);


/**
 * @brief open a file in w write and add mode
 *
 * @param  add [IN]  pathfile
 *
 * @return a FILE*
 */
FILE* openfile_out_add(char* add);


/**
 * @brief open a file in read mode
 *
 * @param  add [IN]  pathfile
 *
 * @return a FILE*
 */
FILE* openfile_in(char* add);

/**
 * @brief print array spin->angles in a FILE
 *
 * @param  spin [IN]  spin grid
 * 
 * @param  fichier [IN]  file where write
 */
void print_spinners(spinners_t* spin, FILE* fichier);


/**
 * @brief print a matrice sizex * sizey, a arrayr of pointeurs
 *
 * @param  matrice [IN] array of pointeur of line of the matrice
 *
 * @param  spin [IN]  spin grid
 * 
 * @param  sizex [IN]  horizontale size of matrice
 * 
 * @param  sizey [IN]  verticale size of matrice
 * 
 * @param  fichier [IN]  file where write
 */
void print_matrice(double** matrice, const int sizex, const int sizey, FILE* fichier);


/**
 * @brief print a matrice_t
 *
 * @param  matrice [IN] ùatrice_t that will be printed
 * 
 * @param  fichier [IN]  file where write
 */
void print_matrice(matrice_t* matrice, FILE* fichier);


/**
 * @brief read a matrice_t
 *
 * @param  matrice [OUT] ùatrice_t that will be read
 * 
 * @param  fichier [OUT]  file where importe
 */
void read_matrice(matrice_t* matrice, char* add);

/**
 * @brief read from a file : array spin->angles and put in a spinners_t, if spin->Ngrid = 1
 *
 * @param  spin [INOUT]  spin grid
 *
 * @param  add [IN]  input file path
 */
void read_spinners(spinners_t* spin, char* add);


/**
 * @brief read from a file : array spin->angles and put in a spinners_t, if spin->Ngrid != 1
 *
 * @param  spin [INOUT]  spin grid
 *
 * @param  add [IN]  input file path
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 * 
 * @param  L [IN]  distance between tow spinners
 */
void read_spinnersall(spinners_t* spin, char* add, int nx, int ny, double L);


/**
 * @brief plot all grid in un spinner_t
 *
 * @param  spin [IN]  spinner_t with Ngrid
 */
void plotall(spinners_t* spin);

/**
 * @brief annealing simualed of a spinner grid
 *
 * @param  spin [INOUT]  spin grid
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 *
 * @param  TO [IN]  initial temparure
 * 
 * @param  TF [IN]  initial temparure
 * 
 * @param  lamnbda [IN]  paramter T *= lambda
 * 
 * @param  Niter [IN]  number of change at a fixed temperature
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 * 
 * @param  seed [IN]  seed for rand_d threat safe
 */
void recuit(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int offset, unsigned int* seed);

/**
 * @brief annealing simualed of a spinner grid
 *
 * @param  spin [INOUT]  spin grid
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 *
 * @param  TO [IN]  initial temparure
 * 
 * @param  TF [IN]  initial temparure
 * 
 * @param  lamnbda [IN]  paramter T *= lambda
 * 
 * @param  Niter [IN]  number of change at a fixed temperature
 * 
 * @param  offest [IN]  position oh the grid considered int spin->angles
 * 
 * @param  seed [IN]  seed for rand_d threat safe
 */
void recuit_anti(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int offset, unsigned int* seed);

/**
 * @brief plot the minimal an maximum energy of a set of spinner lattice
 *
 * @param  add [IN]  add of spinner lattices to be read
 *
 * @param  L [IN]  lenght betwen two spinner center
 *
 * @param  nx [IN]  x-dimension of lattice
 * 
 * @param  ny [IN]  y-dimension of lattice
 */
void plot_MinMax(char* add, double L, int nx, int ny);

/**
 * @brief performe n!/nmin!
 *
 * @param  spin [IN]  spinner_t with Ngrid
 * 
 * @return n!/nmin!
 */
unsigned long factorielle(int n, int nmin);


/**
 * @brief check if 2 grid are equal
 *
 * @param  A [IN]  spinner_t with Ngrid
 * 
 * @param  size [IN]  nx * ny
 * 
 * @param  offset1 [IN]  position of a grid in spin->angles
 *
 * @param  offset2 [IN]  position of a grid in spin->angles
 * 
 * @return returne true if equale, false else
 */
bool isequale(spinners_t* A, const int size, int offset1, int offset2);

/**
 * @brief check if 2 grid are equal
 *
 * @param  A [IN]  spinner_t with Ngrid
 * 
 * @param  size [IN]  nx * ny
 * 
 * @param  offset1 [IN]  position of a grid in spin->angles
 *
 * @param  offset2 [IN]  position of a grid in spin->angles
 * 
 * @return returne true if equale, false else
 */
bool isequale_anti(spinners_t* A, const int size, int offset1, int offset2);


/**
 * @brief remove all duplicate grid of spin->angles
 *
 * @param  spin [INOUT]  spinner_t with Ngrid
 */
void remove_equale(spinners_t* spin);

/**
 * @brief remove all duplicate and no metastable grid of spin->angles
 *
 * @param  spin [INOUT]  spinner_t with Ngrid
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void remove_equale_allmeta(spinners_t* spin, double* H, double* HB);

void remove_equale_allmeta_anti(spinners_t* spin, double* H, double* HB);


/**
 * @brief print the energy of all grid of a spinner_t
 *
 * @param  spin [IN]  spin grid
 *
 * @param  add [IN]  output file path
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void print_E(spinners_t* spin, char* add, double* H, double* HB);


/**
 * @brief plot the mean energy of all grid of a spinner_t and the standard deviation : printf("%f %f %f %f %f\n", track , moyenne, ecart_type, Emin, Emax)
 *
 * @param  spin [IN]  spin grid
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param  track [IN]  double that will plot with the mean energy, that can be use like an abscisse
 */
void plot_E_mean(spinners_t* spin, double* H, double* HB, double track);


/**
 * @brief compute the energy and compute the void print_dist_Histo(spinners_t* spin, char* add, double* H, double* HB, 
					double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB))ograme of all grid of a spinner_t
 *
 * @param  spin [IN]  spin grid
 *
 * @param  add [IN]  output file path
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 */
void print_E_Histo(spinners_t* spin, char* add, double* H, double* HB);

void plot_dist_Histo(spinners_t* spin, FILE* fichier, double* H, double* HB, double track, int Nbin, double min, double max,
					double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB), int p);

void plot_dist_Histo_prob(spinners_t* spin, FILE* fichier, double* H, double* HB, double track, int Nbin, double min, double max,
					double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB), int p);

void plot_stat(spinners_t* spin, double* H, double* HB, double track, 
				double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB) );

void plot_stat_fromone(spinners_t* spin,spinners_t* spin0, double* H, double* HB, double track,
				double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB) );

void to_ppp(spinners_t* spin);

/********************************************************************************/
/*                                                                              */
/*                                distance                                      */
/*                                                                              */
/********************************************************************************/

/**
 * @brief compute the global energy difference :  |E_i-E_j|/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return |E_i-E_j|/N
 */
double dist_EG(spinners_t* spin, int i, int j, int N, double* H, double * HB);

/**
 * @brief compute the local energy difference : \sqrt{\sum_k (E_ik-E_jk)^2}/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return \sqrt{\sum_k (E_ik-E_jk)^2}/N
 */
double dist_EL(spinners_t* spin, int i, int j, int N, double* H, double * HB);

/**
 * @brief compute the local energy difference : \sqrt{\sum_k (E_ik-E_jk)^2}/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return \sqrt{\sum_k (E_ik-E_jk)^2}/N
 */
double dist_EL_anti(spinners_t* spin, int i, int j, int N, double* H, double * HB);

/**
 * @brief compute the Hamming distance between two states : \sqrt{\sum_k (theta_ik-theta_jk)^2}/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return \sqrt{\sum_k (theta_ik-theta_jk)^2}/N
 */
double dist_H(spinners_t* spin, int i, int j, int N, double* H, double * HB);


/**
 * @brief compute a rotation invariant Hamming distance between tow stats : \sqrt{\sum_k \sum_neigbour (theta_ik - theta_ineigbourg[k]-theta_jk + theta_jneigbourg[k])^2}/N
 *
 * @param  spin [IN]  spin grid
 *
 * @param  i [IN]  indice of the i-grid of spin->angles
 * 
 * @param  j [IN]  indice of the j-grid of spin->angles
 * 
 * @param  N [IN]  size of a grid, nx * ny
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @return \sqrt{\sum_k \sum_neigbour (theta_ik - theta_ineigbourg[k]-theta_jk + theta_jneigbourg[k])^2}/N
 */
double dist_HI(spinners_t* spin, int i, int j, int N, double* H, double * HB);

double fromUM(matrice_t* matrice_dist, matrice_t* matrice_ultra);

/**
 * @brief print the matrice of distance bteween two grid of all grid of a spinner_t. The matrice is "trier" and the ultrametric distance
 *
 * @param  spin [IN]  spin grid
 *
 * @param  add [IN]  output file path
 * 
 * @param  distchar [IN]  additional char* that will be insert in add
 * 
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 * 
 * @param dist [in] a fonction that performe the distance betxeen two state
 * 
 * @param  similarity [IN]  the function that compute similarity between cluster
 * 
 * @return the distance betwenn metric and ultrmetric distance
 */
double print_dist(spinners_t* spin, char* add, char* distchar, double*H, double *HB, 
		double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB),
		double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ) );


/********************************************************************************/
/*                                                                              */
/*                                clustering                                    */
/*                                                                              */
/********************************************************************************/

/**
 * @brief compute the distacne between two line of the distance matrice d = sum_l | matrice_il - matrice_jl|
 *
 * @param  B [IN]  a pointeur of the a line of the matrice
 * 
 * @param  A [IN]  a pointeur of the a line of the matrice
 *
 * @param  N [IN]  the size of a line of the matrice 
 * 
 * @return d = sum_l | matrice_il - matrice_jl|
 */
double distline(double* A, double* B, int N);

/**
 * @brief "trie" the distance matrice, by hamming distance on the ultrametric distance
 *
 * @param  matrice [INOUT]  the matrice_t of distance between two states
 * 
 * @param  ultra [INOUT]  the matrice_t of ultrametric distance between two states
 */
void tri(matrice_t* matrice, matrice_t* ultra);

/**
 * @brief compute the average link between 2 cluster  = mean(d(A,B))
 * 
 * @param  cluster_A [IN]  a cluster_t of a tree_t
 * 
 * @param  cluster_B [IN]  a cluster_t of a tree_t
 * 
 * @param  matrice_dist [IN]  the matrice_t of distance of the set of states considered
 * 
 * @return  similarity between two cluster = mean(d(A,B))
 */
double average_link(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist );

/**
 * @brief compute the complete link between 2 cluster  = max(d(A,B))
 * 
 * @param  cluster_A [IN]  a cluster_t of a tree_t
 * 
 * @param  cluster_B [IN]  a cluster_t of a tree_t
 * 
 * @param  matrice_dist [IN]  the matrice_t of distance of the set of states considered
 * 
 * @return  similarity between two cluster = max(d(A,B))
 */
double complete_link(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist );

/**
 * @brief compute the single link between 2 cluster  = min(d(A,B))
 * 
 * @param  cluster_A [IN]  a cluster_t of a tree_t
 * 
 * @param  cluster_B [IN]  a cluster_t of a tree_t
 * 
 * @param  matrice_dist [IN]  the matrice_t of distance of the set of states considered
 * 
 * @return  similarity between two cluster = min(d(A,B))
 */
double single_link(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist );

/**
 * @brief performe the ultrametric distance with an average link clustering
 *
 * @param  tree [IN]  inital partition tree_t of cluster of size 1
 * 
 * @param  matrice_dist [IN]  the matrice_t of distance between two states
 * 
 * @param  matrice_ultra [OUT]  the matrice_t of umtrametric distance between two state
 * 
 * @param  similarity [IN]  the function that compute similarity between cluster
 */
void cluster_fusion(tree_t* tree, matrice_t* matrice_dist, matrice_t* matrice_ultra, double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist )) ;


/**
 * @brief performe the ultrametric distance matrice
 * 
 * @param  matrice_dist [IN]  the matrice_t of distance between two states
 * 
 * @param  matrice_ultra [OUT]  the matrice_t of umtrametric distance between two states
 * 
 * @param  similarity [IN]  the function that compute similarity between cluster
 */
void matrice_ultra(matrice_t* matrice_dist, matrice_t* matrice_ultra, double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ));

/**
 * @brief compute the degree of a spinner lattice
 *
 * @param  spin [IN]  spinner lattice
 * 
 * @param  dist [IN]  a fonction that performe the distance betxeen two state
 * 
 * @param  similarity [IN]  the function that compute similarity between cluster
 * 
 * @return the degree of ultrametricity
 */
long double degree_UM(spinners_t* spin, double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB), 
				double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ) );

void plot_degree_UM_and_stat(spinners_t* spin, double* degree);

void plot_degree_UM_and_stat_anti(spinners_t* spin, double* degree);

/********************************************************************************/
/*                                                                              */
/*                                experimente                                   */
/*                                                                              */
/********************************************************************************/

/**
 * @brief print all metastable state for a n b y lattice
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  L [IN]  lenght betwen two spinner center
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 */
void print_allnxn(char* add, double L, int nx, int ny);

/**
 * @brief print all metastable state for a n x 1 lattice average over disorder between PPM and MPP
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  L [IN]  lenght betwen two spinner center
 * 
 * @param  nl [IN]  n max
 * 
 * @param  Nsimu [IN]  number of simulation for the average
 */
void print_allnx1_recuit(char* add, double L, int nl, int Nsimu);


/**
 * @brief print all metastable state for a 4 by 4 lattice
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  L [IN]  lenght betwen two spinner center
 */
void print_all4x4(char* add, double L);

/**
 * @brief print  metastable state for a n by n lattice
 * 
 * @param  L [IN]  lenght betwen two spinner center
 * 
 * @param  NSimu [IN]  number of simulation for the average
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 */
void print_nxn_anti( double L, int NSimu, int nx, int ny);

/**
 * @brief print all metastable state for a n x 1 lattice average over disorder between PPM and MPP
 * 
 * @param  L [IN]  lenght betwen two spinner center
 * 
 * @param  Nsimu [IN]  number of simulation for the average
 * 
 * @param  nl [IN]  n max
 */
void print_allnx1_anti( double L, int Nsimu, int nl);

/**
 * @brief print all metastable state for a 3 by 3 lattice
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  L [IN]  lenght betwen two spinner center
 */
void print_all3x3(char* add, double L);

/**
 * @brief print all metastable state for a 3 by 3 lattice average over disorder between PPM and MPP
 * 
 * @param  L [IN]  lenght betwen two spinner center
 * 
 * @param  NSimu [IN]  number of simulation for the average
 */
void print_all3x3_anti( double L, int NSimu);

/**
 * @brief print all metastable state for a 2 by 2 lattice
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  L [IN]  lenght betwen two spinner center
 */
void print_all2x2(char* add, double L);

/**
 * @brief print all metastable state for a 2 by 1 lattice average over disorder between PPM and MPP
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  L [IN]  lenght betwen two spinner center
 */
void print_all2x1_anti(char* add, double L);

/**
 * @brief print all metastable state for a 2 by 1 lattice
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  L [IN]  lenght betwen two spinner center
 */
void print_all2x1(char* add, double L);

/**
 * @brief plot all grid of a 3 by 3 lattice
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  L [IN]  lenght betwen two spinner center
 * 
 * @param  bx [IN]  B-field in x direction
 * 
 * @param  by [IN]  B-field in y direction
 */
void plot_all3x3_of_LBX(char* add, double L, double bx , double by, int p );

/**
 * @brief print a the more energtic metastable configuration over Niters rendome itértation
 * 
 * @param  spin [IN]  spin grid
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  Niters [IN]  number of iteration
 */
void print_Emax(spinners_t* spin, char* add, int Niters);

/**
 * @brief print a the less energtic metastable configuration over Niters rendome itértation
 * 
 * @param  spin [IN]  spin grid
 * 
 * @param  add [IN]  pathfile
 * 
 * @param  Niters [IN]  number of iteration
 * 
 * @param  p [IN]  number of threats for OpenMP
 */
void print_Emin( spinners_t* spin,  char* add, int Niters, int p);

/**
 * @brief compute Nsimu annealing of a initial sate
 *
 * @param  spin [INOUT]  spin grid initial and all final grid
 *
 * @param  H [IN]  interaction tensor dipole dipole
 *
 * @param  HB [IN]  interaction tensor dipole B-field
 *
 * @param  TO [IN]  initial temparure
 * 
 * @param  TF [IN]  initial temparure
 * 
 * @param  lamnbda [IN]  paramter T *= lambda
 * 
 * @param  Niter [IN]  number of change at a fixed temperature
 * 
 * @param Nismu [IN] number of annelaing performed
 * 
 * @param p [IN] number of threats for OpenMP
 */
void recuitN(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int Nsimu, int p);

void recuitN_RSBk(spinners_t* spin, double* H, double* HB, int Nrand, double T0, double TF, double T02, double lambda, int Niter, int Nsimu, int NT02);

void recuitN_RSBk_anti(spinners_t* spin, double* H, double* HB, int Nrand, double T0, double TF, double T02, double lambda, int Niter, int Nsimu, int NT02, unsigned int* seed, double prob_anti);

/** 
 * @brief compute degree of ultrametricity of a set of metastable state with average 
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 * 
 * @param  L [IN]  distance between tow spinners
 * 
 * @param  T0 [IN]  first initial temparure
 * 
 * @param  T02 [IN]  second initial temparure
 * 
 * @param  Nrand [IN]  number of random configuration
 * 
 * @param  Niter [IN]  number of ietration for T0
 * 
 * @param  Nmean [IN]  number of configuration for the average
 * 
 * @param  NT02 [IN]  number of iteration for T02
 */
void UM_of_L(int nx, int ny, double L, double T0, double T02, int Nrand, int Niter, int Nmean, int NT02);

/** 
 * @brief compute degree of ultrametricity of a set of metastable state with average and average over quenched random
 * 
 * @param  nx [IN]  x-size of spinner grid
 * 
 * @param  ny [IN]  y-size of spinner grid
 * 
 * @param  L [IN]  distance between tow spinners
 * 
 * @param  T0 [IN]  first initial temparure
 * 
 * @param  T02 [IN]  second initial temparure
 * 
 * @param  Nrand [IN]  number of random configuration
 * 
 * @param  Niter [IN]  number of ietration for T0
 * 
 * @param  Nmean [IN]  number of configuration for the average
 * 
 * @param  NT02 [IN]  number of iteration for T02
 * 
 * @param  prob_anti [IN]  probability of MPP instead of PPM
 * 
 * @param  Nquechedrandom [IN]  number of iteration for the average over quenched random
 */
void UM_of_L_anti(int nx, int ny, double L, double T0, double T02, int Nrand, int Niter, int NT02, double prob_anti, int Nquechedrandom);

/**
 * @brief Recursively generate all configurations with a given number of -1s in a 9-element array.
 *
 * @param index Start index for placing -1s.
 * @param remaining Number of -1s to place.
 * @param configurations Output array of configurations.
 * @param count Pointer to number of configurations generated.
 * @param a Temporary array used to build configurations.
 */
void insert_ones(int index, int remaining, int configurations[512][9], int *count, int a[9]);

/**
 * @brief Generate all 9-element configurations with exactly n values set to -1.
 *
 * @param n Number of -1s.
 * @param configurations Output array of configurations.
 * @param count Pointer to the number of configurations generated.
 */
void generate_configurations(int n, int configurations[512][9], int *count);

/**
 * @brief Print all generated configurations.
 *
 * @param configurations Configuration array.
 * @param count Number of configurations.
 */
void print_configurations(int configurations[512][9], int count);

/**
 * @brief Evaluate all configurations of 9 spinners and compute statistics with average over queched disorder.
 *
 * @param L Lattice parameter.
 */
void anti9_all(double L);


