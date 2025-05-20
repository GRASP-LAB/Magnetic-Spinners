#include "discret_model_function.h"
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <vector>
#include <omp.h>
#include <climits>
#include <float.h>


/********************************************************************************/
/*                                                                              */
/*                            comput Energy dip�le fct                          */
/*                                                                              */
/********************************************************************************/

double Udipole(double m1_x, double m1_y, double m2_x, double m2_y, double r1_x, double r1_y, double r2_x, double r2_y) {

	double dx = r1_x - r2_x; 
	double dy = r1_y - r2_y; 
	double r = sqrt(pow(dx, 2) + pow(dy, 2)); 
	return MU0 * 0.25 / PI * (m1_x * m2_x + m1_y * m2_y - 3. * (m1_x * dx + m1_y * dy) * (m2_x * dx + m2_y * dy) / pow(r, 2))
		/ pow(r, 3) ;
}

double UB(double mx, double my, double bx, double by) {
	return (-mx * bx - my * by) / HREF / MU;
}

double UBspinner(int theta, double bx, double by) {
	double U = 0.;
	for (int i = 0; i < 6; i +=2) {
		double ang = PI3 * (theta + i);
		double Udd = UB(cos(ang), sin(ang), bx, by);
		if (i == 4) { Udd *= -1.; }
		U += Udd;
	}
	return U;
}

double Uspinner(int theta, int thetav, double L, int nv) {

	double U = 0.;
	for (int i = 0; i < 6; i += 2) {
		for (int j = 0; j < 6; j += 2) {

			double angv = PI3 * (thetav + j); 
			double ang = PI3 * (theta + i); 

			double Udd = Udipole(cos(ang), sin(ang), cos(angv), sin(angv),
				R * cos(ang), R * sin(ang), L * cos(PI3 * nv) + R * cos(angv), L * sin(PI3 * nv) + R * sin(angv));

			if ((j == 4 && i != 4) || (j != 4 && i == 4)) { Udd *= -1.; } // The minus sign is the dipole at j = 4, which is the bottom left
			U += Udd; 
		}
	}
	return U / HREF;
}

double Uspinner_pertubed(int theta, int thetav, double L, int nv, double beta) {

	double U = 0.;
	for (int i = 0; i < 6; i += 2) {
		for (int j = 0; j < 6; j += 2) {

			double angv = PI3 * (thetav + j); 
			double ang = PI3 * (theta + i) + beta * PI / 180.; 

			double Udd = Udipole(cos(ang), sin(ang), cos(angv), sin(angv),
				R * cos(ang), R * sin(ang), L * cos(PI3 * nv) + R * cos(angv), L * sin(PI3 * nv) + R * sin(angv));

			if ((j == 4 && i != 4) || (j != 4 && i == 4)) { Udd *= -1.; } // The minus sign is the dipole at j = 4, which is the bottom left
			U += Udd; 
		}
	}
	return U / HREF;
}

double Uspinner_ppp(int theta, int thetav, double L, int nv) {

	double U = 0.;
	for (int i = 0; i < 6; i += 2) {
		for (int j = 0; j < 6; j += 2) {

			double angv = PI3 * (thetav + j); 
			double ang = PI3 * (theta + i); 

			double Udd = Udipole(cos(ang), sin(ang), cos(angv), sin(angv),
				R * cos(ang), R * sin(ang), L * cos(PI3 * nv) + R * cos(angv), L * sin(PI3 * nv) + R * sin(angv));
			U += Udd; 
		}
	}
	return U / HREF;
}


double*  H_init(double L) {

	double* H = (double*)malloc(SIZE_H * sizeof(double));

	if (H == NULL) {
		fprintf(stderr, "H_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 6; i++) {
		for (int v = 0; v < 6; v++) {
			SETELEMENT(H, i, v, Uspinner(i, v, L, 0));
		}
	}
	return H;
}

double*  H_init_L(double L) {

	double* H = (double*)malloc(SIZE_H * sizeof(double));

	if (H == NULL) {
		fprintf(stderr, "H_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}

	double ref = abs(Udipole(1,0,1,0,R,0,L-R,0));

	for (int i = 0; i < 6; i++) {
		for (int v = 0; v < 6; v++) {
			SETELEMENT(H, i, v, Uspinner(i, v, L, 0) / ref * HREF);
		}
	}

	return H;
}


double*  H_init_L25(double L) {

	double* H = (double*)malloc(SIZE_H * sizeof(double));

	if (H == NULL) {
		fprintf(stderr, "H_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}

	double ref = abs(Udipole(1,0,1,0,R,0,L-R,0));

	double MAX = - DBL_MAX;

	for (int i = 0; i < 6; i++) {
		for (int v = 0; v < 6; v++) {
			SETELEMENT(H, i, v, Uspinner(i, v, L, 0) / ref * HREF);
			if(MAX < abs(GETELEMENT(H,i,v)) ) MAX = abs(GETELEMENT(H,i,v));
		}
	}

	for (int i = 0; i < 6; i++) {
		for (int v = 0; v < 6; v++) {
			SETELEMENT(H, i, v, GETELEMENT(H,i,v) / MAX * H25MAX );
		}
	}

	return H;
}

double*  H_ppp_init(double L) {

	double* H = (double*)malloc(SIZE_H * sizeof(double));

	if (H == NULL) {
		fprintf(stderr, "H_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 6; i++) {
		for (int v = 0; v < 6; v++) {
			SETELEMENT(H, i, v, Uspinner_ppp(i, v, L, 0));
		}
	}
	return H;
}

double* H_B_init( double bx, double by) {

	double* HB = (double*)malloc(6 * sizeof(double));

	if (HB == NULL) {
		fprintf(stderr, "H_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}


	for (int i = 0; i < 6; i++) { HB[i] = UBspinner(i, bx, by); }
	return HB;
}

double* H_B_init_L( double bx, double by, double L) {

	double* HB = (double*)malloc(6 * sizeof(double));
	double ref = Udipole(1,0,1,0,R,0,L-R,0);
	if (HB == NULL) {
		fprintf(stderr, "H_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}


	for (int i = 0; i < 6; i++) { HB[i] = UBspinner(i, bx, by) / ref * HREF; }
	return HB;
}

void H_plot(double *H) {

	for (int i = 0; i < 6; i++) {
		printf("%f", GETELEMENT(H, i, 0));
		for (int v = 1; v < 6; v++) {
			printf("\t %f", GETELEMENT(H, i, v));
		}
		printf("\n");
	}
	printf("\n\n");
}

void H_B_plot(double* HB) {

	printf("%f", HB[0]);
	for (int i = 1; i < 6; i++) {
		printf("\t %f", HB[i]);
	}
	printf("\n\n");
}



/********************************************************************************/
/*                                                                              */
/*          Initialisation spinner lattice and lattice fct                      */
/*                                                                              */
/********************************************************************************/

void spinners_init(spinners_t *spin, double L, int nx, int ny, int Ngrid){

	spin->L = L;
	spin->nx = nx;
	spin->ny = ny;
	spin->Ngrid = Ngrid;

	spin->angles = (int*)malloc(nx * ny * Ngrid * sizeof(int));

	if (spin->angles == NULL) {
		fprintf(stderr, "spinner_init : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}
}

void spinners_init_anti(spinners_t *spin, double L, int nx, int ny, unsigned int* seed, int Nanti){

	spin->L = L;
	spin->nx = nx;
	spin->ny = ny;
	spin->Ngrid = 1;

	spin->angles = (int*)malloc(nx * ny * sizeof(int));
	spin->anti = (int*)malloc(nx * ny * sizeof(int));

	if (spin->angles == NULL) {
		fprintf(stderr, "spinner_init_anti : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}
	if (spin->anti == NULL) {
		fprintf(stderr, "spinner_init_anti : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}

	if(Nanti == 0) for(int j = 0; j < nx * ny; j++) spin->anti[j] = 1 ;
	if(Nanti == nx * ny) for(int j = 0; j < nx * ny; j++) spin->anti[j] = -1 ;
	else{
		for(int j = 0; j < nx * ny; j++) spin->anti[j] = 1 ;
		int Nrand = 0;
		int* rand = (int*)malloc(sizeof(int) * nx * ny);
		for(int i = 0; i < nx * ny; i++) rand[i] = 100;
		rand[0] = rand_r(seed) % (nx * ny);
		while(Nrand!= Nanti){
			int index = rand_r(seed) % (nx * ny);
			bool alerady = false;
			for(int i = 0; i <= Nrand ; i++){
				if(rand[i] == index ) alerady = true;
			}
			if(!alerady){
				Nrand++;
				rand[Nrand] = index;
			}
		
		}
		for(int j = 0; j < Nanti; j++) spin->anti[rand[j]] = -1;
		free(rand);
	}
}

void spinners_init_anti_prob(spinners_t *spin, double L, int nx, int ny, int Ngrid, unsigned int* seed, double prob_anti){

	spin->L = L;
	spin->nx = nx;
	spin->ny = ny;
	spin->Ngrid = 1;

	spin->angles = (int*)malloc(nx * ny * sizeof(int) * Ngrid);
	spin->anti = (int*)malloc(nx * ny * sizeof(int) * Ngrid);

	if (spin->angles == NULL) {
		fprintf(stderr, "spinner_init_anti : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}
	if (spin->anti == NULL) {
		fprintf(stderr, "spinner_init_anti : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < spin->nx * spin->ny * Ngrid; i++) {
		if ((double)rand_r(seed) / (double)RAND_MAX < prob_anti ) spin->anti[i] = -1;
		else spin->anti[i] = 1;
	}

}

void Finalisation_simu(spinners_t* spin, double *H, double* HB) {
	free(spin->angles);
	free(H);
	free(HB);
}

int* neighbour( spinners_t* spin, int index) { //  Starting from the top left corner with kx = ky = 0

    int ky = (index - index % spin->nx) / spin->nx;
	int kx = index - ky * spin->nx; 

	int* voisins = (int*)malloc(SIZE_NEIGHBOUR * sizeof(int));

	if (voisins == NULL) {
		fprintf(stderr, "neighbour : Allocation de m�moire �chou�e.\n");
		exit(EXIT_FAILURE);
	}


    if (kx != 0) {
        voisins[3] = ky * spin->nx + kx - 1; // left
    }
	else { voisins[3]  = - 1; }

    if (ky != 0)     // top
    {
        if (ky % 2 == 1)     // top left
        {
            voisins[2] = (ky - 1) * spin->nx + kx;
        }
        else if (ky % 2 == 0 && kx > 0) {
            voisins[2] = (ky - 1) * spin->nx + kx - 1;
        }
		else {
			voisins[2] = -1;
		}

        if (ky % 2 == 0)     // top right
        {
            voisins[1] = (ky - 1) * spin->nx + kx;
        }
        else if (ky % 2 == 1 && kx < spin->nx - 1) {
            voisins[1] = (ky - 1) * spin->nx + kx + 1;
		}
		else { voisins[1] = -1; }
    }
	else {
		voisins[2] = -1;
		voisins[1] = -1;
	}

    if (kx != spin->nx - 1) {
        voisins[0] = ky * spin->nx + kx + 1;    //right
	}
	else { voisins[0] = -1; }

    if (ky != spin->ny - 1)     //bottom
    {
        if (ky % 2 == 1)     // bottom left
        {
			voisins[4] = (ky + 1) * spin->nx + kx; 
        }
        else if (ky % 2 == 0 && kx > 0) {
            voisins[4] = (ky + 1) * spin->nx + kx - 1; 
		}
		else { voisins[4] = -1; }

        if (ky % 2 == 0)     // bottom right
        {
            voisins[5] = (ky + 1) * spin->nx + kx;
        }
        else if (ky % 2 == 1 && kx < spin->nx - 1) {
            voisins[5] = (ky + 1) * spin->nx + kx + 1;
		}
		else { voisins[5] = -1; }
	} else {// Order: 0 right, 1 top-right, 2 top-left, 3 left, 4 bottom-left, 5 bottom-right
		voisins[4] = -1;
		voisins[5] = -1;
	}
	return voisins;
}

/********************************************************************************/
/*                                                                              */
/*                           magnetisation                                      */
/*                                                                              */
/********************************************************************************/

double magnetisation(spinners_t* spin){
	double mx = 0;
	double my = 0;
	int N = spin->nx * spin->ny;
	for(int i = 0; i < N; i++){
		mx += cos(PI3 * (spin->angles[i] - 2));
		my += sin(PI3 * (spin->angles[i] - 2));
	}

	return sqrt(mx * mx + my * my) / (double)N;
}


/********************************************************************************/
/*                                                                              */
/*                           Energy of spinner lattice                          */
/*                                                                              */
/********************************************************************************/

double E_local(spinners_t* spin, int index, double* H, double* HB, int offset) {

	int* voisins = neighbour(spin, index);
	double U = 0.;
	for (int v = 0; v < SIZE_NEIGHBOUR; v++) {
		if (voisins[v] != -1) {
			U += GETELEMENT(H, (spin->angles[index + offset] - v + 6) % 6, (spin->angles[voisins[v] + offset] - v + 6) % 6); 
		}
	}
	free(voisins); 
	return U + HB[( spin->angles[index] + 6) % 6];
}

double E_local_anti(spinners_t* spin, int index, double* H, double* HB, int offset) {

	int* voisins = neighbour(spin, index);
	double U = 0.;
	for (int v = 0; v < SIZE_NEIGHBOUR; v++) {
		if (voisins[v] != -1) {
			U += GETELEMENT(H, (spin->angles[index + offset] - v + 6) % 6, (spin->angles[voisins[v] + offset] - v + 6) % 6) * spin->anti[index + offset] * spin->anti[voisins[v] + offset];
		}
	}
	free(voisins);
	return U + HB[( spin->angles[index] + 6) % 6];
}

double E_total(spinners_t* spin, double*H, double* HB, int offset) {
	double U = 0.;
	for (int j = 0; j < spin->ny; j++) 
	{
		for (int i = 0; i < spin->nx; i ++) 
		{
			U += E_local(spin, i + j * spin->nx, H, HB, offset);
		}
	}
	return U / 2.;
}

double E_total_anti(spinners_t* spin, double*H, double* HB, int offset) {
	double U = 0.;
	for (int j = 0; j < spin->ny; j++) 
	{	
		for (int i = 0; i < spin->nx; i ++) 
		{
			U += E_local_anti(spin, i + j * spin->nx, H, HB, offset);
			
		}
	}
	return U / 2.;
}

bool metastable(spinners_t* spin, double* H, double* HB, int offset) {
	int N = spin->ny * spin->nx;
	for(int i = 0; i < N; i++){
		spin->angles[i]++;
		double EP = E_local(spin, i, H, HB, offset);
		spin->angles[i] -= 2;
		double EM = E_local(spin, i, H, HB, offset);
		spin->angles[i]++;
		double ER = E_local(spin, i, H, HB, offset);
		if (EM + DBL_PRECISION < ER || EP + DBL_PRECISION < ER) { 
			return false;
		}
	}
	return true;
}

bool metastable_anti(spinners_t* spin, double* H, double* HB, int offset) {
	int N = spin->ny * spin->nx;
	for(int i = 0; i < N; i++){
		spin->angles[i]++;
		double EP = E_local_anti(spin, i, H, HB, offset);
		spin->angles[i] -= 2;
		double EM = E_local_anti(spin, i, H, HB, offset);
		spin->angles[i]++;
		double ER = E_local_anti(spin, i, H, HB, offset);
		if (EM + DBL_PRECISION < ER || EP + DBL_PRECISION < ER) { 
			return false;
		}
	}
	return true;
}

/********************************************************************************/
/*                                                                              */
/*         						 to continu     			                    */
/*                                                                              */
/********************************************************************************/

bool meta_continu(spinners_t* spin, double* H, double* HB,  int offset){
	int N = spin->nx * spin->ny;
	for(int i = 0; i < N; i++){
		int* voisins = neighbour(spin, i);
		double ER = E_local(spin, i, H, HB, offset);
		double EP = 0;
		double EM = 0;
		for(int j = 0; j < SIZE_NEIGHBOUR; j++){
			if(voisins[j] != -1){
				EP += Uspinner_pertubed(spin->angles[i + offset], spin->angles[voisins[j] + offset], spin->L, j, META_ANGLE);
				EM += Uspinner_pertubed(spin->angles[i + offset], spin->angles[voisins[j] + offset], spin->L, j, - META_ANGLE);
			}
		}
		
		if (EM + DBL_PRECISION < ER || EP + DBL_PRECISION < ER) { 
			return false;
		}
		free(voisins);
	}
	return true;
}

void to_meta_continu(spinners_t* spin){
	double* H = H_init_L25(spin->L);
	double* HB = H_B_init_L(0, 0, spin->L);
	int N = spin->nx * spin->ny;

	bool* meta = (bool*)malloc(spin->Ngrid * sizeof(bool));
	int Nmeta = 0;

	for(int i = 0; i < spin->Ngrid; i++){
		if(meta_continu(spin, H, HB, i * N)){
			meta[i] = true;
			Nmeta++;
		}
		else meta[i] = false;
	}

	int* spin_buffer = (int*)malloc(Nmeta * N * sizeof(int));

	int k = 0;
	for(int i = 0; i < spin->Ngrid; i++){
		if(meta[i]){
			for(int j = 0; j < N; j++){
				spin_buffer[k * N + j] = spin->angles[ i * N + j];
			}
			k++;
		}
	}

	int* buffer = spin_buffer;
	spin_buffer = spin->angles;
	spin->angles = spin_buffer;
	spin->Ngrid = Nmeta;

	free(meta);
	free(spin_buffer);
	free(H);
	free(HB);
}


/********************************************************************************/
/*                                                                              */
/*                                      utilitaire                              */
/*                                                                              */
/********************************************************************************/

void change(int* A, int* B){
	int c = *A;
	*A = *B;
	*B = c;
}

double mean(double* A, int N){
	double mean = 0;
	for(int i = 0; i < N; i++) mean += A[i];
	return mean / (double)N;
}

double SD(double* A, int N){
	double meanA = mean(A, N);
	double sd = 0;
	for(int i = 0; i < N; i++) sd += (A[i] - meanA) * (A[i] - meanA);
	return sd / (double)(N - 1);
}

FILE* openfile_out(char* add) {

	 FILE* fichier = fopen(add, "w");

    if (fichier == NULL) {
        fprintf(stderr, "openfile_out : Impossible d'ouvrir le fichier %s pour l'�criture.\n", add);
        exit(EXIT_FAILURE);
    }
	return fichier;
}

FILE* openfile_out_add(char* add) {

	 FILE* fichier = fopen(add, "a");

    if (fichier == NULL) {
        fprintf(stderr, "openfile_out : Impossible d'ouvrir le fichier %s pour l'�criture.\n", add);
        exit(EXIT_FAILURE);
    }
	return fichier;
}

FILE* openfile_in(char* add) {

	FILE* fichier = fopen( add, "r");

	if (fichier == NULL) {
		fprintf(stderr, "openfile_in : Impossible d'ouvrir le fichier %s pour lecture.\n", add);
		exit(EXIT_FAILURE);
	}
	return fichier;
}


void print_spinners(spinners_t* spin, FILE* fichier) {

	int N = spin->nx * spin->ny;
	for(int j = 0; j < spin->Ngrid; j++){
		fprintf(fichier, "%d", spin->angles[0]);
		for (int i = 1; i < N; i++) {
			fprintf(fichier, "\t%d", (spin->angles[i + N * j] + 6) % 6);
		}
		fprintf(fichier, "\n");
	}
}

void print_matrice(double** matrice, const int sizex, const int sizey, FILE* fichier) {

	for(int i = 0; i < sizey ; i++){
		fprintf(fichier, "%f", matrice[i][0]);
		for(int j = 1; j < sizex ; j++){
			fprintf(fichier, "\t%f", matrice[i][j]);
		}
		if(i < sizey - 1){fprintf(fichier, "\n");}
	}
}


void print_matrice(matrice_t* matrice, FILE* fichier) {

	for(int i = 0; i < matrice->N ; i++){
		fprintf(fichier, "%f", matrice->line[i].col[0]);
		for(int j = 1; j < matrice->N ; j++){
			fprintf(fichier, "\t%f", matrice->line[i].col[j]);
		}
		if(i < matrice->N - 1){fprintf(fichier, "\n");}
	}
}

void read_matrice(matrice_t* matrice, char* add){
	FILE* fichier = openfile_in(add);

	std::vector<double> tab(0);
	double value;
	while (fscanf(fichier, "%lf", &value) == 1) {tab.push_back(value);}

	int N = sqrt(tab.size());
	if(N * N != tab.size()) printf("read_matrice, size did not matche %s \n", add);

	matrice->N = N;
	matrice->line = (matrice_line_t*)malloc(N * sizeof(matrice_line_t));
	if (matrice->line == NULL) {
		fprintf(stderr, "read_matrice, matrice : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < N; i++){
		matrice->line[i].pos = i;
		matrice->line[i].col = (double*)malloc(N * sizeof(double));
		if (matrice->line[i].col == NULL) {
		fprintf(stderr, "read_matrice, matrice[%d].col : Allocation de memoire echouee.\n", i);
		exit(EXIT_FAILURE);
		}
	}

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			matrice->line[j].col[i] = tab[i * N + j];
		}
	}

	fclose(fichier);
}

void read_spinners(spinners_t* spin, char* add) {

	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}
	FILE* fichier = openfile_in(add);
	int N = spin->nx * spin->ny;
	for (int i = 0; i < N; i++) {fscanf(fichier, "%d", &spin->angles[i]); }
	fclose(fichier);
}

void read_spinnersall(spinners_t* spin, char* add, int nx, int ny, double L) {

	FILE* fichier = openfile_in(add);
	int N = nx * ny;
	int value;
	std::vector<int> input(0);
	while (fscanf(fichier, "%d", &value) == 1) {input.push_back(value);}
	spin->Ngrid = (int)input.size() / nx / ny; 
	free(spin->angles);
	spin->angles = (int*)malloc(input.size() * sizeof(int));
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < N; j++){
			spin->angles[N * i + j] = input[N * i + j];
		}
	}
	fclose(fichier);
}

void plotall(spinners_t* spin){
	int N = spin->nx * spin->ny;
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < N; j++){
			printf("%d\t", spin->angles[i * N + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void recuit(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int offset, unsigned int* seed) {
	int N = spin->ny * spin->nx;
	for (double T = T0; T >= TF; T *= lambda) {
		for (int i = 0; i < Niter; i++) {
			int index = rand_r(seed) % N;
			double E0 = E_local(spin, index, H, HB, offset);

			int signe = rand_r(seed) % 2 == 0 ? 1 : -1;
			spin->angles[index + offset] += signe;  
			double EF = E_local(spin, index, H, HB, offset); 
			if (E0 + DBL_PRECISION < EF) { 
				double temp = (double)rand_r(seed) / (double)RAND_MAX;
				if (temp > exp((E0 - EF) / T)) { spin->angles[index + offset] -= signe; }	
			}
			spin->angles[index + offset] = (spin->angles[index + offset] + 6) % 6;
		}
	}
}

void recuit_anti(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int offset, unsigned int* seed) {
	int N = spin->ny * spin->nx;
	for (double T = T0; T >= TF; T *= lambda) {
		for (int i = 0; i < Niter; i++) {
			int index = rand_r(seed) % N;
			double E0 = E_local_anti(spin, index, H, HB, offset);

			int signe = rand_r(seed) % 2 == 0 ? 1 : -1;
			spin->angles[index + offset] += signe;  
			double EF = E_local_anti(spin, index, H, HB, offset); 
			if (E0 + DBL_PRECISION < EF) { 
				double temp = (double)rand_r(seed) / (double)RAND_MAX;
				if (temp > exp((E0 - EF) / T)) { spin->angles[index + offset] -= signe; }	
			}
			spin->angles[index + offset] = (spin->angles[index + offset] + 6) % 6;
		}
	}
}

void plot_MinMax(char* add, double L, int nx, int ny){
	spinners_t spin;
	spinners_init(&spin,  L, nx, ny, 1);
	spinners_t spinmin;
	spinners_init(&spinmin,  L, nx, ny, 1);
	spinners_t spinmax;
	spinners_init(&spinmax,  L, nx, ny, 1);
	double* H = H_init(L);
	double* HB = H_B_init(0, 0);
	const int N = nx * ny;

	read_spinnersall(&spin, add, nx, ny, L);

	double Emax = -DBL_MAX;
	double Emin = DBL_MAX;

	for(int i = 0; i < spin.Ngrid; i++){
		double E = E_total(&spin, H, HB, N * i);
		if(E < Emin){
			Emin = E;
			for(int j = 0; j < N; j++){spinmin.angles[j] = spin.angles[N * i + j];}
		}
		if(E > Emax){
			Emax = E;
			for(int j = 0; j < N; j++){spinmax.angles[j] = spin.angles[N * i + j];}
		}
	}
	printf("Emin %f\n\n", Emin);
	for(int j = 0; j < N; j++){printf("%d\t", spinmin.angles[j]);}
	printf("\n\nEmax %f\n\n", Emax);
	for(int j = 0; j < N; j++){printf("%d\t", spinmax.angles[j]);}
	printf("\n\n");
	Finalisation_simu(&spin, H, HB);
	free(spinmax.angles);
	free(spinmin.angles);
}

unsigned long factorielle(int n, int nmin) {
	if (n == 0 || n == 1) {return 1; }
	if (nmin == n) { return nmin; }
	else {
		return n * factorielle(n - 1, nmin);
	}
}

bool isequale(spinners_t* A, const int size, int offset1, int offset2){ 
	for(int i = 0; i < size; i++){
		if (A->angles[offset1 + i] != A->angles[offset2 + i]){
			return false;
		}
	}
	return true;
}

bool isequale_anti(spinners_t* A, const int size, int offset1, int offset2){ 
	for(int i = 0; i < size; i++){
		if (A->angles[offset1 + i] != A->angles[offset2 + i] || A->anti[offset1 + i] != A->anti[offset2 + i]){
			return false;
		}
	}
	return true;
}

void remove_equale(spinners_t* spin){
	int N = spin->nx * spin->ny;
	bool* unicity = (bool*)malloc(spin->Ngrid * sizeof(bool));
	if (unicity == NULL) {
		fprintf(stderr, "remove_equal, unicity : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < spin->Ngrid; i++){unicity[i] = true;}
	int sizeout = spin->Ngrid;
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int  j = i + 1; j < spin->Ngrid; j++){ 
				if(isequale(spin, N, i * N,  j * N) && unicity[i]){
					unicity[i] = false;
					sizeout--;
				}
			}
		}
	}

	int* out = (int*)malloc(sizeout * N * sizeof(int));
	if (out == NULL) {
		fprintf(stderr, "remove_equal, out : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	int k = 0;
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int j = 0; j < N ; j++){
				out[N * k + j] = spin->angles[N * i + j];
			}
			k++;
		}
	}

	spin->angles = (int*)realloc(spin->angles, sizeout * N * sizeof(int));
	if (spin->angles == NULL) {
        fprintf(stderr, "remove_equal, new_angles :Réallocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	
	memcpy(spin->angles, out, sizeout * N * sizeof(int));
	free(out);
	free(unicity);
	printf("Ngird = %d, after remove_equale %d\n", spin->Ngrid, sizeout );
	spin->Ngrid = sizeout;
}

void remove_equale_allmeta(spinners_t* spin, double* H, double* HB){
	int N = spin->nx * spin->ny;
	bool* unicity = (bool*)malloc(spin->Ngrid * sizeof(bool));
	if (unicity == NULL) {
		fprintf(stderr, "remove_equal_allmeta, unicity : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < spin->Ngrid; i++){unicity[i] = true;}
	int sizeout = spin->Ngrid;
	for(int i =0; i < spin->Ngrid; i++){
		if(!metastable(spin, H, HB, N*i)){
			unicity[i] = false;
			sizeout--;
		};
	}
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int  j = i + 1; j < spin->Ngrid; j++){ 
				if(isequale(spin, N, i * N,  j * N) && unicity[i]){
					unicity[i] = false;
					sizeout--;
				}
			}
		}
	}

	int* out = (int*)malloc(sizeout * N * sizeof(int));
	if (out == NULL) {
		fprintf(stderr, "remove_equal_allmeta, out : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	int k = 0;
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int j = 0; j < N ; j++){
				out[N * k + j] = spin->angles[N * i + j];
			}
			k++;
		}
	}

	spin->angles = (int*)realloc(spin->angles, sizeout * N * sizeof(int));
	if (spin->angles == NULL) {
        fprintf(stderr, "remove_equal_allmeta, new_angles :Réallocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	
	memcpy(spin->angles, out, sizeout * N * sizeof(int));
	free(out);
	free(unicity);
	
	spin->Ngrid = sizeout;
}

void remove_equale_allmeta_anti(spinners_t* spin, double* H, double* HB){
	int N = spin->nx * spin->ny;
	bool* unicity = (bool*)malloc(spin->Ngrid * sizeof(bool));
	if (unicity == NULL) {
		fprintf(stderr, "remove_equal_allmeta, unicity : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < spin->Ngrid; i++){unicity[i] = true;}
	int sizeout = spin->Ngrid;
	for(int i =0; i < spin->Ngrid; i++){
		if(!metastable_anti(spin, H, HB, N*i)){
			unicity[i] = false;
			sizeout--;
		};
	}
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int  j = i + 1; j < spin->Ngrid; j++){ 
				if(isequale_anti(spin, N, i * N,  j * N) && unicity[i]){
					unicity[i] = false;
					sizeout--;
				}
			}
		}
	}

	int* out = (int*)malloc(sizeout * N * sizeof(int));
	if (out == NULL) {
		fprintf(stderr, "remove_equal_allmeta, out : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	int* out_anti = (int*)malloc(sizeout * N * sizeof(int));
	if (out_anti == NULL) {
		fprintf(stderr, "remove_equal_allmeta, out : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	int k = 0;
	for(int i = 0; i < spin->Ngrid; i++){
		if(unicity[i]){
			for(int j = 0; j < N ; j++){
				out[N * k + j] = spin->angles[N * i + j];
				out_anti[N * k + j] = spin->anti[N * i + j];
			}
			k++;
		}
	}

	spin->angles = (int*)realloc(spin->angles, sizeout * N * sizeof(int));
	spin->anti = (int*)realloc(spin->anti, sizeout * N * sizeof(int));
	if (spin->angles == NULL || spin->anti == NULL) {
        fprintf(stderr, "remove_equal_allmeta, new_angles :Réallocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	
	memcpy(spin->angles, out, sizeout * N * sizeof(int));
	free(out);
	free(unicity);
	
	spin->Ngrid = sizeout;
}

void print_E(spinners_t* spin, char* add, double* H, double* HB){
	FILE* fichier = openfile_out(add);
	int N = spin->nx * spin->ny;
	double* E = (double*)malloc(spin->Ngrid * sizeof(double));
	if (E == NULL) {
        fprintf(stderr, "print_E, E : allocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	for(int i = 0; i < spin->Ngrid; i++){ E[i] = E_total(spin, H, HB, N * i); }
	for(int i = 0; i < spin->Ngrid - 1; i++){ fprintf(fichier, "%f\n", E[i]); }
	fprintf(fichier, "%f", E[spin->Ngrid - 1]);
	free(E);
	fclose(fichier);
}

void plot_E_mean(spinners_t* spin, double* H, double* HB, double track){
	int N = spin->nx * spin->ny;
	double* E = (double*)malloc(spin->Ngrid * sizeof(double));
	if (E == NULL) {
        fprintf(stderr, "print_E, E : allocation de mémoire échouée.\n");
        exit(EXIT_FAILURE);
    }
	double Emin = DBL_MAX;
	double Emax = -DBL_MAX;
	for(int i = 0; i < spin->Ngrid; i++){ 
		E[i] = E_total(spin, H, HB, N * i);
		if(E[i] > Emax) Emax = E[i];
		if(E[i] < Emin) Emin = E[i];
	}
	double moyenne = 0.0;
    for (int i = 0; i < spin->Ngrid; i++) { moyenne += E[i]; }
    moyenne /= spin->Ngrid;
    double SD = 0.0;
    for (int i = 0; i < spin->Ngrid; i++) { SD += (E[i] - moyenne)*(E[i] - moyenne);}
    SD = sqrt(SD / spin->Ngrid);

	double skewness = 0;
	double kurtosis = 0;

	for(int i = 0; i < spin->Ngrid ; i++ ){
		double tmp = (moyenne - E[i]);
		skewness += tmp * tmp * tmp;
		kurtosis += tmp * tmp * tmp * tmp;
	}

	skewness /= (double)(N - 1) * SD * SD * SD;
	kurtosis /= (double)(N - 1) * SD * SD * SD * SD;

	printf("%f %f %f %f %f %f %f\n", track , moyenne, SD, Emin, Emax, skewness, kurtosis);
	free(E);
}


void print_E_Histo(spinners_t* spin, char* add, double* H, double* HB){
	FILE* fichier = openfile_out(add);
	int N = spin->nx * spin->ny;
	std::map<double, int> histogram;
	for(int i = 0; i < spin->Ngrid; i++){ histogram[E_total(spin, H, HB, N * i)]++; }
	for (const auto& entry : histogram){ fprintf(fichier, "%f %d\n", entry.first, entry.second); }
	fclose(fichier);
}

void plot_dist_Histo(spinners_t* spin, FILE* fichier, double* H, double* HB, double track, int Nbin, double min, double max,
					double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB), int p){
	
	int N = spin->nx * spin->ny;
	double Lbin = (max - min) / (double)Nbin;
	int* histo = (int*)calloc(Nbin, sizeof(int));

	#pragma omp parallel for num_threads(p) collapse(2)
	for(int i = 0; i < spin->Ngrid; i++) {
    	for(int j = i + 1; j < spin->Ngrid; j++) {
        	double distance = dist(spin, i, j, N, H, HB);
        	int index = (int)((distance - min) / Lbin);
        	
			#pragma omp atomic
            histo[index]++;
    	}
	}

	fprintf(fichier, "%f", track);
	for(int i = Nbin - 1; i >= 0; i--){ fprintf(fichier, "\t%d", histo[i]); }
	fprintf(fichier, "\n");
	fflush(fichier);
	
	free(histo);
}

void plot_dist_Histo_prob(spinners_t* spin, FILE* fichier, double* H, double* HB, double track, int Nbin, double min, double max,
					double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB), int p){
	
	int N = spin->nx * spin->ny;
	double Lbin = (max - min) / (double)Nbin;
	double* histo = (double*)calloc(Nbin, sizeof(double));

	#pragma omp parallel for num_threads(p)  collapse(2)
	for(int i = 0; i < spin->Ngrid; i++) {
    	for(int j = i + 1; j < spin->Ngrid; j++) {
        	double distance = dist(spin, i, j, N, H, HB);
        	int index = (int)((distance - min) / Lbin);
        	
			#pragma omp atomic
            histo[index]+=1.;
    	}
	}

	double norm = 0;
	for(int i = 0; i < Nbin; i++)  norm += histo[i];
	for(int i = 0; i < Nbin; i++)  histo[i] /= norm;

	fprintf(fichier, "%f", track);
	for(int i = Nbin - 1; i >= 0; i--){ fprintf(fichier, "\t%f", histo[i]); }
	fprintf(fichier, "\n");
	fflush(fichier);
	
	free(histo);
}

void plot_stat(spinners_t* spin, double* H, double* HB, double track,
				double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB) ){
	
	
	int N = spin->Ngrid * (spin->Ngrid + 1) / 2 - spin->Ngrid;
	int n = spin->nx * spin->ny;
	int k = 0;

	double min = DBL_MAX;
	double max = -DBL_MAX;

	double* matrice = (double*)malloc( N * sizeof(double));
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < i; j++){
			matrice[k] = dist( spin, i, j, n, H, HB);
			if(matrice[k] < min) min = matrice[k];
			if(matrice[k] > max) max = matrice[k] ;
			k++;
		}
	}
	
	remove_equale_allmeta(spin, H, HB);

	int Nmeta = spin->Ngrid;

	double mean = 0;
	double SD = 0;
	double skewness = 0;
	double kurtosis = 0;

	for(int i = 0; i < N ; i++ ){
		mean += matrice[i];
	}
	mean /= (double)N;
	for(int i = 0; i < N ; i++ ){
		double tmp = (mean - matrice[i]);
		SD += tmp * tmp;
		skewness += tmp * tmp * tmp;
		kurtosis += tmp * tmp * tmp * tmp;
	}
	SD /= (double)(N - 1);
	SD = sqrt(SD);

	skewness /= (double)(N - 1) * SD * SD * SD;
	kurtosis /= (double)(N - 1) * SD * SD * SD * SD;

	printf("%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", track, Nmeta, min, max, mean, SD, skewness, kurtosis);
	free(matrice);
}

void plot_stat_fromone(spinners_t* spin,spinners_t* spin0, double* H, double* HB, double track,
				double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB) ){
	
	remove_equale_allmeta(spin, H, HB);

	int Nmeta = spin->Ngrid;

	
	int N = spin->Ngrid ;
	int n = spin->nx * spin->ny;

	double min = DBL_MAX;
	double max = -DBL_MAX;

	spin->angles = (int*)realloc(spin->angles, sizeof(int) * n * (N + 1));

	for(int i = 0; i < n; i++) {spin->angles[ spin->Ngrid * n + i] = spin0->angles[i];}
	spin->Ngrid++;

	double* matrice = (double*)malloc( N * sizeof(double));
	for(int i = 0; i < N; i++){
		matrice[i] = dist( spin, i, N, n, H, HB);
		if(matrice[i] < min) min = matrice[i];
		if(matrice[i] > max) max = matrice[i] ;
	}

	double mean = 0;
	double SD = 0;
	double skewness = 0;
	double kurtosis = 0;

	for(int i = 0; i < N ; i++ ){
		mean += matrice[i];
	}
	mean /= (double)N;
	for(int i = 0; i < N ; i++ ){
		double tmp = (mean - matrice[i]);
		SD += tmp * tmp;
		skewness += tmp * tmp * tmp;
		kurtosis += tmp * tmp * tmp * tmp;
	}
	SD /= (double)(N - 1);
	SD = sqrt(SD);

	skewness /= (double)(N - 1) * SD * SD * SD;
	kurtosis /= (double)(N - 1) * SD * SD * SD * SD;

	printf("%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", track, Nmeta, min, max, mean, SD, skewness, kurtosis);
	free(matrice);
}


void to_ppp(spinners_t* spin){
	const int N = spin->nx * spin->ny;
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < N; j++){
			spin->angles[i * N + j] = ((spin->angles[i * N + j] + 6 ) % 6) % 2 ;
		}
	}
}

/********************************************************************************/
/*                                                                              */
/*                                distance                                      */
/*                                                                              */
/********************************************************************************/

double dist_EG(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	return abs(E_total(spin, H, HB, N * i) - E_total(spin, H, HB, N * j)) / (double)N;
}

double dist_EL(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	double d = 0.;
	for(int k = 0; k < N; k++){
    double E = E_local(spin, k, H, HB, i * N) - E_local(spin, k, H, HB, j * N);
    d += E *E ;
  }
	return sqrt(d / (double)N);
}

double dist_EL_anti(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	double d = 0.;
	for(int k = 0; k < N; k++){
    double E = E_local_anti(spin, k, H, HB, i * N) - E_local_anti(spin, k, H, HB, j * N);
    d += E *E ;
  }
	return sqrt(d / (double)N);
}

double dist_H(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	double d = 0.;
	for(int k = 0; k < N; k++){
   int delta = (spin->angles[i * N + k] - spin->angles[j * N + k]) * (spin->angles[i * N + k] - spin->angles[j * N + k]);
    switch( delta ){
      case 25:
        d += 1.;
        break;
      case 16 :
        d += 4.;
        break;
      default:
        d += delta;
        break;
    }
  }
	return sqrt(d / (double)N) ;
}

double dist_HI(spinners_t* spin, int i, int j, int N, double* H, double * HB){
	double d = 0.;
	for(int k = 0; k < N; k++){
		double dd = 0;
		int* voisins = neighbour(spin, k);
		for(int l = 0; l < SIZE_NEIGHBOUR; l++){
			if(voisins[l] != -1){
        
       int delta = spin->angles[i * N + k] - spin->angles[i * N + voisins[l]] - spin->angles[j * N + k] + spin->angles[j * N + voisins[l]];
       delta *= delta;
       switch( delta ){
          case 25:
          dd += 1.;
          break;
        case 16 :
          dd += 4.;
          break;
        default:
          dd += delta;
          break;
        }
      }
		}
		d += dd;
		free(voisins);
	}
	return sqrt(d / (double)N) ;
}

double fromUM(matrice_t* matrice_dist, matrice_t* matrice_ultra){
	double d = 0;
	for(int i = 0; i < matrice_ultra->N ; i++){
		for(int j = i + 1; j < matrice_ultra->N ; j++){
			d += (matrice_dist->line[i].col[j] - matrice_ultra->line[i].col[j]) * (matrice_dist->line[i].col[j] - matrice_ultra->line[i].col[j]);
		}
	}
	return sqrt(d / (double)matrice_ultra->N);
}

double print_dist(spinners_t* spin, char* add, char* distchar, double*H, double *HB,
		double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB), 
		double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist )){
	char path[STRING_MAX];
	strcpy(path, add);
	strcat(path, "_L"); 
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", spin->L);
	strcat(path, "_dist"); 
	strcat(path, distchar); 
	strcat(path, ".txt");


	char path2[STRING_MAX];
	strcpy(path2, add);
	strcat(path2, "_L"); 
	snprintf(path2 + strlen(path2), sizeof(path2) - strlen(path2), "%f", spin->L);
	strcat(path2, "_dist"); 
	strcat(path2, distchar); 
	strcat(path2, "_ultra.txt");

	FILE* fichier = openfile_out(path);
	FILE* fichier_ultra = openfile_out(path2);

	matrice_t matrice;
	matrice_t ultra;
	matrice.N = spin->Ngrid;
	ultra.N = spin->Ngrid;
	matrice.line = (matrice_line_t*)malloc( matrice.N * sizeof(matrice_line_t));
	ultra.line = (matrice_line_t*)malloc( ultra.N * sizeof(matrice_line_t));
	if (matrice.line == NULL || ultra.line == NULL) {
		fprintf(stderr, "print_matrice, matrice.line and ultra.line : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < matrice.N; i++){
		matrice.line[i].col = (double*)malloc( matrice.N * sizeof(double));
		ultra.line[i].col = (double*)malloc( ultra.N * sizeof(double));
		if (matrice.line[i].col == NULL || ultra.line[i].col == NULL) {
			fprintf(stderr, "print_matrice, matrice.line[%d].col and ultra.line[%d].col : Allocation de memoire echouee.\n", i, i);
			exit(EXIT_FAILURE);
		}
		matrice.line[i].pos = i;
		ultra.line[i].pos = i;
	}

	const int N = spin->nx * spin->ny;
	for(int i = 0; i < spin->Ngrid ; i ++){
		for(int j = i; j < spin->Ngrid ; j ++){
			double d = dist(spin, i, j, N, H, HB);
			matrice.line[i].col[j] = d;
			matrice.line[j].col[i] = d;
		}
	}
	
	matrice_ultra(&matrice, &ultra, similarity);
	tri(&matrice, &ultra);
	print_matrice(&matrice, fichier);
	print_matrice(&ultra, fichier_ultra);
	double deltaUM =fromUM(&matrice, &ultra);
	fclose(fichier);
	fclose(fichier_ultra);
	for(int i = 0; i < matrice.N; i++){free(matrice.line[i].col);}
	free(matrice.line);
	for(int i = 0; i < ultra.N; i++){free(ultra.line[i].col);}
	free(ultra.line);
	return deltaUM;
}

/********************************************************************************/
/*                                                                              */
/*                                clustering                                    */
/*                                                                              */
/********************************************************************************/

double distline(double* A, double* B, int N){
	double d = 0;
	for(int i = 0; i < N; i++){ d += fabs(A[i] - B[i])*fabs(A[i] - B[i]);}
	return sqrt(d);
}

void tri(matrice_t* matrice, matrice_t* ultra){
	
	for (int i = 0; i < matrice->N - 1; i++) {
        double dmin = distline(matrice->line[i].col, matrice->line[i + 1].col, matrice->N);
        for (int j = i + 2; j < matrice->N; j++) {
			double d = distline(matrice->line[i].col, matrice->line[j].col, matrice->N);
            if (d < dmin ) {
				double* change = matrice->line[j].col;
				matrice->line[j].col = matrice->line[i + 1].col;
				matrice->line[i + 1].col = change; 

				dmin = d;

				int temp = matrice->line[i + 1].pos;
				matrice->line[i + 1].pos = matrice->line[j].pos;
				matrice->line[j].pos = temp;
            }
        }
    }

	
	double* change1 = (double*)malloc(matrice->N * sizeof(double));
	double* change2 = (double*)malloc(matrice->N * sizeof(double));
	double** change3 = (double**)malloc(matrice->N * sizeof(double*));
	if (change1 == NULL || change2 == NULL || change3 == NULL) {
		fprintf(stderr, "tri, change : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < matrice->N; i++){change3[i] = ultra->line[i].col;}
	for(int i = 0; i < matrice->N; i++){ultra->line[i].col = change3[matrice->line[i].pos];}
	free(change3);

	for(int i = 0; i < ultra->N; i++){
		memcpy(change1, matrice->line[i].col, matrice->N * sizeof(double));
		memcpy(change2, ultra->line[i].col, ultra->N * sizeof(double));
		for(int j = 0; j < ultra->N; j++){ultra->line[i].col[j] = change2[matrice->line[j].pos];}
		for(int j = 0; j < matrice->N; j++){matrice->line[i].col[j] = change1[matrice->line[j].pos];}
	}
	free(change1);
	free(change2);

}

double average_link(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ){
	double similarity = 0.;
	for(int k = 0; k < cluster_A.size;  k++){
		for(int l = 0; l < cluster_B.size;  l++){
			similarity += matrice_dist->line[cluster_A.pos[k]].col[cluster_B.pos[l]];
		}
	}
	return similarity / (cluster_A.size * cluster_B.size);
}

double complete_link(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ){
	double similarity = -DBL_MAX;
	for(int k = 0; k < cluster_A.size;  k++){
		for(int l = 0; l < cluster_B.size;  l++){
			double d = matrice_dist->line[cluster_A.pos[k]].col[cluster_B.pos[l]];
			if(similarity < d) similarity = d;
		}
	}
	return similarity ;
}

double single_link(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ){
	double similarity = DBL_MAX;
	for(int k = 0; k < cluster_A.size;  k++){
		for(int l = 0; l < cluster_B.size;  l++){
			double d = matrice_dist->line[cluster_A.pos[k]].col[cluster_B.pos[l]];
			if(similarity > d) similarity = d;
		}
	}
	return similarity;
}


void cluster_fusion(tree_t* tree, matrice_t* matrice_dist, matrice_t* matrice_ultra, double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist )) {
	
	int Ncluster = tree->N;
	int* newClusterPos = (int*)malloc(tree->N * sizeof(int));
    if (newClusterPos == NULL) {
        fprintf(stderr, "cluster_fusion, newClusterPos : Allocation de memoire echouee.\n");
        exit(EXIT_FAILURE);
    }

    while (Ncluster > 1) {
        
        int minI = -1, minJ = -1;
        double minSimilarite = DBL_MAX;
		double similarity_dbl = 0.;

        for (int i = 0; i < tree->N; i++) {
            if (tree->cluster[i].notmerged) {
                for (int j = i + 1; j < tree->N; j++) {
                    if (tree->cluster[j].notmerged) {
                        similarity_dbl = similarity(tree->cluster[i],tree->cluster[j], matrice_dist);
                        if ( similarity_dbl < minSimilarite) {
                            minSimilarite = similarity_dbl;
                            minI = i;
                            minJ = j;
                        }
                    }
                }
            }
        }

		for (int k = 0; k < tree->cluster[minI].size; k++) {
            for (int l = 0; l < tree->cluster[minJ].size; l++) {
				matrice_ultra->line[tree->cluster[minI].pos[k]].col[tree->cluster[minJ].pos[l]] = minSimilarite;
				matrice_ultra->line[tree->cluster[minJ].pos[l]].col[tree->cluster[minI].pos[k]] = minSimilarite;
        	}
        }

        int newClusterSize = tree->cluster[minI].size + tree->cluster[minJ].size;

        memcpy(newClusterPos, tree->cluster[minI].pos, tree->cluster[minI].size * sizeof(int));
        memcpy(newClusterPos + tree->cluster[minI].size, tree->cluster[minJ].pos, tree->cluster[minJ].size * sizeof(int));

        tree->cluster[minI].size = newClusterSize;
        tree->cluster[minI].pos = (int*)realloc(tree->cluster[minI].pos, newClusterSize * sizeof(int));
		if (tree->cluster[minI].pos == NULL) {
            fprintf(stderr, "cluster_fusion, tree->cluster[%d].pos : Allocation de memoire echouee.\n", minI);
            exit(EXIT_FAILURE);
        }
		memcpy(tree->cluster[minI].pos, newClusterPos, newClusterSize * sizeof(int));

        tree->cluster[minJ].notmerged = false;

		Ncluster--;
    }
	free(newClusterPos);
}

void matrice_ultra(matrice_t* matrice_dist, matrice_t* matrice_ultra, double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist )){
	tree_t tree;
	tree.N = matrice_dist->N;
	tree.cluster = (cluster_t*)malloc( tree.N * sizeof(cluster_t));
	if (tree.cluster == NULL) {
		fprintf(stderr, "matrice_ultra, tree.cluster : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < tree.N; i++){
		tree.cluster[i].pos = (int*)malloc(sizeof(int));
		if (tree.cluster[i].pos == NULL) {
			fprintf(stderr, "matrice_ultra, tree.cluster[%d].pos : Allocation de memoire echouee.\n", i);
			exit(EXIT_FAILURE);
		}
		tree.cluster[i].size = 1;
		tree.cluster[i].pos[0] = i;
		tree.cluster[i].notmerged = true;
	}

	for(int i = 0; i < matrice_ultra->N; i++){matrice_ultra->line[i].col[i] = 0;}

	cluster_fusion(&tree, matrice_dist, matrice_ultra, similarity);

	for(int i = 0; i < tree.N; i++){free(tree.cluster[i].pos);}
	free(tree.cluster);
}

long double degree_UM(spinners_t* spin, double (*dist)(spinners_t* spin, int i, int j, int N, double* H, double * HB), 
				double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ) ){
	long double sum_diff = 0;
	long double sum_metric = 0;
	const int N = spin->nx * spin->ny;
	double* H = H_init(spin->L);
	double* HB = H_B_init(0, 0);

	matrice_t matrice;
	matrice_t ultra;
	matrice.N = spin->Ngrid;
	ultra.N = spin->Ngrid;
	matrice.line = (matrice_line_t*)malloc( matrice.N * sizeof(matrice_line_t));
	ultra.line = (matrice_line_t*)malloc( ultra.N * sizeof(matrice_line_t));
	if (matrice.line == NULL || ultra.line == NULL) {
		fprintf(stderr, "print_matrice, matrice.line and ultra.line : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < matrice.N; i++){
		matrice.line[i].col = (double*)malloc( matrice.N * sizeof(double));
		ultra.line[i].col = (double*)malloc( ultra.N * sizeof(double));
		if (matrice.line[i].col == NULL || ultra.line[i].col == NULL) {
			fprintf(stderr, "print_matrice, matrice.line[%d].col and ultra.line[%d].col : Allocation de memoire echouee.\n", i, i);
			exit(EXIT_FAILURE);
		}
		matrice.line[i].pos = i;
		ultra.line[i].pos = i;
	}

	for(int i = 0; i < spin->Ngrid ; i ++){
		for(int j = i; j < spin->Ngrid ; j ++){
			double d = dist(spin, i, j, N, H, HB);
			matrice.line[i].col[j] = d;
			matrice.line[j].col[i] = d;
		}
	}

	matrice_ultra(&matrice, &ultra, similarity);
	
	for(int i = 0; i < spin->Ngrid ; i ++){
		for(int j = 0; j < i ; j ++){
			sum_diff += matrice.line[i].col[j] - ultra.line[i].col[j];
			sum_metric += matrice.line[i].col[j];
		}
	}

	free(H);
	free(HB);
	for(int i = 0; i < matrice.N; i++){free(matrice.line[i].col);}
	free(matrice.line);
	for(int i = 0; i < ultra.N; i++){free(ultra.line[i].col);}
	free(ultra.line);
	return abs(sum_diff)/ sum_metric;
}

void plot_degree_UM_and_stat(spinners_t* spin, double* degree){
	
	const int N = spin->nx * spin->ny;

	matrice_t matriceEG;
	matrice_t ultraEG;
	matriceEG.N = spin->Ngrid;
	ultraEG.N = spin->Ngrid;
	matrice_t matriceEL;
	matrice_t ultraEL;
	matriceEL.N = spin->Ngrid;
	ultraEL.N = spin->Ngrid;
	matrice_t matriceH;
	matrice_t ultraH;
	matriceH.N = spin->Ngrid;
	ultraH.N = spin->Ngrid;
	matrice_t matriceIH;
	matrice_t ultraIH;
	matriceIH.N = spin->Ngrid;
	ultraIH.N = spin->Ngrid;
	matriceEG.line = (matrice_line_t*)malloc( matriceEG.N * sizeof(matrice_line_t));
	ultraEG.line = (matrice_line_t*)malloc( ultraEG.N * sizeof(matrice_line_t));
	matriceEL.line = (matrice_line_t*)malloc( matriceEG.N * sizeof(matrice_line_t));
	ultraEL.line = (matrice_line_t*)malloc( ultraEG.N * sizeof(matrice_line_t));
	matriceH.line = (matrice_line_t*)malloc( matriceEG.N * sizeof(matrice_line_t));
	ultraH.line = (matrice_line_t*)malloc( ultraEG.N * sizeof(matrice_line_t));
	matriceIH.line = (matrice_line_t*)malloc( matriceEG.N * sizeof(matrice_line_t));
	ultraIH.line = (matrice_line_t*)malloc( ultraEG.N * sizeof(matrice_line_t));
	if (matriceEG.line == NULL || ultraEG.line == NULL || matriceEL.line == NULL || ultraEL.line == NULL 
	|| matriceH.line == NULL || ultraH.line == NULL || matriceIH.line == NULL || ultraIH.line == NULL) {
		fprintf(stderr, "print_matrice, matrice.line and ultra.line : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < matriceEG.N; i++){
		matriceEG.line[i].col = (double*)malloc( matriceEG.N * sizeof(double));
		ultraEG.line[i].col = (double*)malloc( ultraEG.N * sizeof(double));
		matriceEL.line[i].col = (double*)malloc( matriceEG.N * sizeof(double));
		ultraEL.line[i].col = (double*)malloc( ultraEG.N * sizeof(double));
		matriceH.line[i].col = (double*)malloc( matriceEG.N * sizeof(double));
		ultraH.line[i].col = (double*)malloc( ultraEG.N * sizeof(double));
		matriceIH.line[i].col = (double*)malloc( matriceEG.N * sizeof(double));
		ultraIH.line[i].col = (double*)malloc( ultraEG.N * sizeof(double));
		if (matriceEG.line[i].col == NULL || ultraEG.line[i].col == NULL || matriceEL.line[i].col == NULL || ultraEL.line[i].col == NULL 
		|| matriceH.line[i].col == NULL || ultraH.line[i].col == NULL || matriceIH.line[i].col == NULL || ultraIH.line[i].col == NULL) {
			fprintf(stderr, "print_matrice, matrice.line[%d].col and ultra.line[%d].col : Allocation de memoire echouee.\n", i, i);
			exit(EXIT_FAILURE);
		}
		matriceEG.line[i].pos = i;
		ultraEG.line[i].pos = i;
		matriceEL.line[i].pos = i;
		ultraEL.line[i].pos = i;
		matriceH.line[i].pos = i;
		ultraH.line[i].pos = i;
		matriceIH.line[i].pos = i;
		ultraIH.line[i].pos = i;
	}

	
	#pragma omp parallel for
	for(int i = 0; i < spin->Ngrid ; i ++){

		double* H = H_init_L25(spin->L);
		double* HB = H_B_init(0, 0);

		for(int j = i; j < spin->Ngrid ; j ++){
			double dEG = dist_EG(spin, i, j, N, H, HB);
			double dEL = dist_EL(spin, i, j, N, H, HB);
			double dH = dist_H(spin, i, j, N, H, HB);
			double dIH = dist_HI(spin, i, j, N, H, HB);


			matriceEG.line[i].col[j] = dEG;
			matriceEG.line[j].col[i] = dEG;
			matriceEL.line[i].col[j] = dEL;
			matriceEL.line[j].col[i] = dEL;
			matriceH.line[i].col[j] = dH;
			matriceH.line[j].col[i] = dH;
			matriceIH.line[i].col[j] = dIH;
			matriceIH.line[j].col[i] = dIH;
		}

		free(H);
		free(HB);
	}

	if(omp_get_num_threads() >= 4){
		#pragma omp parallel
		{
			switch (omp_get_thread_num()) {
            case 0:
                matrice_ultra(&matriceEG, &ultraEG, single_link);
                break;
            case 1:
                matrice_ultra(&matriceEL, &ultraEL, single_link);
                break;
            case 2:
                matrice_ultra(&matriceH, &ultraH, single_link);
                break;
            case 3:
                matrice_ultra(&matriceIH, &ultraIH, single_link);
                break;
            default:
                break;
        	}
		}
	} 
	else {
		matrice_ultra(&matriceEG, &ultraEG, single_link);
		matrice_ultra(&matriceEL, &ultraEL, single_link);
		matrice_ultra(&matriceH, &ultraH, single_link);
		matrice_ultra(&matriceIH, &ultraIH, single_link);
	}

	tri(&ultraEL, &matriceEL);
	
	FILE* file = openfile_out("/mnt/c/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/recuit_of_L/test.txt");
	print_matrice(&matriceEL, file);
	fclose(file);

	FILE* fileU = openfile_out("/mnt/c/Users/axelf/OneDrive - Universite de Liege/Mémoire/simulation/recuit_of_L/testU.txt");
	print_matrice(&ultraEL, fileU);
	fclose(fileU);

	long double sum_diffEG = 0;
	long double sum_metricEG = 0;
	long double sum_diffEL = 0;
	long double sum_metricEL = 0;
	long double sum_diffH = 0;
	long double sum_metricH = 0;
	long double sum_diffIH = 0;
	long double sum_metricIH = 0;

	#pragma omp parallel for collapse(2) reduction(+:sum_diffEG) reduction(+:sum_metricEG)\
										 reduction(+:sum_diffEL) reduction(+:sum_metricEL)\
										 reduction(+:sum_diffH) reduction(+:sum_metricH)\
										 reduction(+:sum_diffIH) reduction(+:sum_metricIH)
	for(int i = 0; i < matriceEG.N ; i ++){
		for(int j = 0; j < matriceEG.N ; j ++){
			sum_diffEG += matriceEG.line[i].col[j] - ultraEG.line[i].col[j];
			sum_metricEG += matriceEG.line[i].col[j];
			sum_diffEL += matriceEL.line[i].col[j] - ultraEL.line[i].col[j];
			sum_metricEL += matriceEL.line[i].col[j];
			sum_diffH += matriceH.line[i].col[j] - ultraH.line[i].col[j];
			sum_metricH += matriceH.line[i].col[j];
			sum_diffIH += matriceIH.line[i].col[j] - ultraIH.line[i].col[j];
			sum_metricIH += matriceIH.line[i].col[j];
		}
	}

	for(int i = 0; i < matriceEG.N; i++){
		free(matriceEG.line[i].col);
		free(matriceEL.line[i].col);
		free(matriceH.line[i].col);
		free(matriceIH.line[i].col);
		free(ultraEG.line[i].col);
		free(ultraEL.line[i].col);
		free(ultraH.line[i].col);
		free(ultraIH.line[i].col);
	}
	free(matriceEG.line);
	free(matriceEL.line);
	free(matriceH.line);
	free(matriceIH.line);
	free(ultraEG.line);
	free(ultraEL.line);
	free(ultraH.line);
	free(ultraIH.line);


	degree[0] = abs(sum_diffEG)/ sum_metricEG;
	degree[1] = abs(sum_diffEL)/ sum_metricEL;
	degree[2] = abs(sum_diffH)/ sum_metricH;
	degree[3] = abs(sum_diffIH)/ sum_metricIH;
}

void plot_degree_UM_and_stat_anti(spinners_t* spin, double* degree){
	
	const int N = spin->nx * spin->ny;

	matrice_t matriceEL;
	matrice_t ultraEL;
	matriceEL.N = spin->Ngrid;
	ultraEL.N = spin->Ngrid;
	
	matriceEL.line = (matrice_line_t*)malloc( matriceEL.N * sizeof(matrice_line_t));
	ultraEL.line = (matrice_line_t*)malloc( ultraEL.N * sizeof(matrice_line_t));

	if (matriceEL.line == NULL || ultraEL.line == NULL) {
		fprintf(stderr, "print_matrice, matrice.line and ultra.line : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < matriceEL.N; i++){
		matriceEL.line[i].col = (double*)malloc( matriceEL.N * sizeof(double));
		ultraEL.line[i].col = (double*)malloc( ultraEL.N * sizeof(double));
		if (matriceEL.line[i].col == NULL || ultraEL.line[i].col == NULL) {
			fprintf(stderr, "print_matrice, matrice.line[%d].col and ultra.line[%d].col : Allocation de memoire echouee.\n", i, i);
			exit(EXIT_FAILURE);
		}
		matriceEL.line[i].pos = i;
		ultraEL.line[i].pos = i;
	}

	
	#pragma omp parallel for
	for(int i = 0; i < spin->Ngrid ; i ++){

		double* H = H_init_L25(spin->L);
		double* HB = H_B_init(0, 0);

		for(int j = i; j < spin->Ngrid ; j ++){
			double dEL = dist_EL_anti(spin, i, j, N, H, HB);

			matriceEL.line[i].col[j] = dEL;
			matriceEL.line[j].col[i] = dEL;
		}

		free(H);
		free(HB);
	}

	matrice_ultra(&matriceEL, &ultraEL, single_link);

	tri(&ultraEL, &matriceEL);


	long double sum_diffEL = 0;
	long double sum_metricEL = 0;

	#pragma omp parallel for collapse(2) reduction(+:sum_diffEL) reduction(+:sum_metricEL)
	for(int i = 0; i < matriceEL.N ; i ++){
		for(int j = 0; j < matriceEL.N ; j ++){
			
			sum_diffEL += matriceEL.line[i].col[j] - ultraEL.line[i].col[j];
			sum_metricEL += matriceEL.line[i].col[j];
		}
	}

	for(int i = 0; i < matriceEL.N; i++){
		free(matriceEL.line[i].col);
		free(ultraEL.line[i].col);
	}

	free(ultraEL.line);
	free(matriceEL.line);


	*degree = abs(sum_diffEL)/ sum_metricEL;
}
/********************************************************************************/
/*                                                                              */
/*                                experimente                                   */
/*                                                                              */
/********************************************************************************/

void print_allnxn(char* add, double L, int nx, int ny){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;

	if(nx * ny > 16) printf("lattice too large\n");

	int size[16];
	for(int i = 0; i < 16; i++){
		if(i + 1 <= nx * ny) size[i] = 6;
		else size[i] = 1;
	}

  #pragma omp parallel num_threads(1)
  {
  	spinners spin;
	spinners_init(&spin, L, 16, 1, 1);
	spin.nx = nx;
	spin.ny = ny;
	#pragma omp for collapse(6)
	for(int i0 = 0 ; i0 < size[0]; i0++){
		for(int i1 = 0 ; i1 < size[1]; i1++){
			for(int i2 = 0 ; i2 < size[2]; i2++){
				for(int i3 = 0 ; i3 < size[3]; i3++){
					for(int i4 = 0 ; i4 < size[4]; i4++){
						for(int i5 = 0 ; i5 < size[5]; i5++){
							spin.angles[0] = i0;
							spin.angles[1] = i1;
							spin.angles[2] = i2;
							spin.angles[3] = i3;
							spin.angles[4] = i4;
							spin.angles[5] = i5;
							for(int i6 = 0 ; i6 < size[6]; i6++){
								spin.angles[6] = i6;
								for(int i7 = 0 ; i7 < size[7]; i7++){
									spin.angles[7] = i7;
									for(int i8 = 0 ; i8 <size[8]; i8++){
										spin.angles[8] = i8;
										for(int i9 = 0 ; i9 < size[9]; i9++){
											spin.angles[9] = i9;
											for(int i10 = 0 ; i10 < size[10]; i10++){
												spin.angles[10] = i10;
												for(int i11 = 0 ; i11 < size[11]; i11++){
													spin.angles[11] = i11;
													for(int i12 = 0 ; i12 < size[12]; i12++){
														spin.angles[12] = i12;
														for(int i13 = 0 ; i13 < size[13]; i13++){
															spin.angles[13] = i13;
															for(int i14 = 0 ; i14 < size[14]; i14++){
																spin.angles[14] = i14;
																for(int i15 = 0 ; i15 < size[15]; i15++){
																	spin.angles[15] = i15;
																
																	if(metastable(&spin, H, HB, 0)) {

																		print_spinners(&spin, fichier);
																		
																		#pragma omp atomic   
																		n++;
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
  free(spin.angles);
  }
	
	printf("%f\t%d\t%d\t%d\n",L, nx, ny, n);
	free(H);
	free(HB);
	fclose(fichier);
}

void print_allnx1_anti( double L, int Nsimu, int nl){
	
	double* H = H_init_L(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;

	int size[16];
	for(int i = 0; i < 16; i++){
		if(i + 1 <= nl) size[i] = 6;
		else size[i] = 1;
	
	}

	double* E = (double*)malloc(sizeof(double) * Nsimu * (nl + 1));
	int* N = (int*)malloc(sizeof(int) * Nsimu * (nl + 1));

  #pragma omp parallel num_threads(10)
  {

	unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();

	#pragma omp for collapse(2)
	for(int l = 0; l < nl + 1; l++){
	for(int ns = 0; ns < Nsimu + 1; ns++){
	int n = 0;
  	spinners spin;
	spinners_init_anti(&spin, L, nl, 1, &seed, l);
	int* buff = (int*)malloc(sizeof(int) * 16);
	int* temp = spin.angles;
	for(int h = 0; h < nl; h++) buff[h] = spin.angles[h];
	spin.angles = buff; 
	spin.nx = nl;
	double e = 0;
	for(int i0 = 0 ; i0 < size[0]; i0++){
		for(int i1 = 0 ; i1 < size[1]; i1++){
			for(int i2 = 0 ; i2 < size[2]; i2++){
				for(int i3 = 0 ; i3 < size[3]; i3++){
					for(int i4 = 0 ; i4 < size[4]; i4++){
						for(int i5 = 0 ; i5 < size[5]; i5++){
							spin.angles[0] = i0;
							spin.angles[1] = i1;
							spin.angles[2] = i2;
							spin.angles[3] = i3;
							spin.angles[4] = i4;
							spin.angles[5] = i5;
							for(int i6 = 0 ; i6 < size[6]; i6++){
								spin.angles[6] = i6;
								for(int i7 = 0 ; i7 < size[7]; i7++){
									spin.angles[7] = i7;
									for(int i8 = 0 ; i8 <size[8]; i8++){
										spin.angles[8] = i8;
										for(int i9 = 0 ; i9 < size[9]; i9++){
											spin.angles[9] = i9;
											for(int i10 = 0 ; i10 < size[10]; i10++){
												spin.angles[10] = i10;
												for(int i11 = 0 ; i11 < size[11]; i11++){
													spin.angles[11] = i11;
													for(int i12 = 0 ; i12 < size[12]; i12++){
														spin.angles[12] = i12;
														for(int i13 = 0 ; i13 < size[13]; i13++){
															spin.angles[13] = i13;
															for(int i14 = 0 ; i14 < size[14]; i14++){
																spin.angles[14] = i14;
																for(int i15 = 0 ; i15 < size[15]; i15++){
																	spin.angles[15] = i15;
																	spin.nx = nl;
																	if(metastable_anti(&spin, H, HB, 0)) {

																		double e_temp = E_total_anti(&spin, H, HB, 0);
																		#pragma omp atomic
																		e += e_temp;
																		
																		n++;
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	E[l * Nsimu + ns] = e /(double)n / (double)nl;
	N[l * Nsimu + ns] = n;
  free(spin.angles);
  free(spin.anti);
  free(temp);

	}
  }
  }

	for(int i =0; i < nl + 1; i++){
	
	double meanE =0;
	double SDE = 0;
	double meanN =0;
	double SDN = 0;

	for(int j = 0; j < Nsimu; j++){
		meanE += E[i * Nsimu + j];
		meanN += N[i * Nsimu + j];
	}

	meanE /= (double)Nsimu;
	meanN /= (double)Nsimu;

	for(int j = 0; j < Nsimu; j++){
		SDE += (meanE - E[i * Nsimu + j]) * (meanE - E[i * Nsimu + j]);
		SDN += (meanN - N[i * Nsimu + j]) * (meanN - N[i * Nsimu + j]);
	}

  	printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n", nl, 1, Nsimu, i, meanE, sqrt(SDE / (double)Nsimu), meanN, sqrt(SDN / (double)Nsimu));
	}

	free(E);
	free(N);
	free(H);
	free(HB);
	
}

void print_allnx1_recuit(char* add, double L, int nl, int Nsimu){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	
	spinners spin;
	spinners_init(&spin, L, nl, 1, Nsimu);


	srand(time(NULL));
	unsigned int seed = time(NULL) + 1000;

	for(int i = 0; i < Nsimu; i++){
		for(int j = 0; j < nl; j++) spin.angles[nl * i + j] = rand() % 6;	
		recuit(&spin, H, HB, 0.000001, 0.000000001, 0.75, nl * 5, i, &seed);												
	}

	remove_equale_allmeta(&spin, H, HB);

	double* E = (double*)malloc(Nsimu * sizeof(double));

	double mean = 0;
	double SD = 0;

	for(int i = 0; i < spin.Ngrid; i++){
		double e = E_total(&spin, H, HB, nl * i) / (double)nl;
		fprintf(fichier, "%d\t%f\n", nl, e);
		E[i] = e;
		mean += e;
	}

	for(int i = 0; i < spin.Ngrid; i++){
		SD += (mean - E[i]) * (mean - E[i]);
	}
	
	SD = sqrt(SD / (double)spin.Ngrid);

	free(E);
	free(spin.angles);
	printf("%d\t%d\t%f\t%f\n", nl, spin.Ngrid, mean / (double)spin.Ngrid, SD );
	free(H);
	free(HB);
	fclose(fichier);
}

void print_all4x4(char* add, double L){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;

	
  #pragma omp parallel
  {
  	spinners spin;
	spinners_init(&spin, L, 4, 4, 1);
	#pragma omp for collapse(6)
	for(int i0 = 0 ; i0 < 6; i0++){
		for(int i1 = 0 ; i1 < 6; i1++){
			for(int i2 = 0 ; i2 < 6; i2++){
				for(int i3 = 0 ; i3 < 6; i3++){
					for(int i4 = 0 ; i4 < 6; i4++){
						for(int i5 = 0 ; i5 < 6; i5++){
							spin.angles[0] = i0;
							spin.angles[1] = i1;
							spin.angles[2] = i2;
							spin.angles[3] = i3;
							spin.angles[4] = i4;
							spin.angles[5] = i5;
							for(int i6 = 0 ; i6 < 6; i6++){
								spin.angles[6] = i6;
								for(int i7 = 0 ; i7 < 6; i7++){
									spin.angles[7] = i7;
									for(int i8 = 0 ; i8 < 6; i8++){
										spin.angles[8] = i8;
										for(int i9 = 0 ; i9 < 6; i9++){
											spin.angles[9] = i9;
											for(int i10 = 0 ; i10 < 6; i10++){
												spin.angles[10] = i10;
												for(int i11 = 0 ; i11 < 6; i11++){
													spin.angles[11] = i11;
													for(int i12 = 0 ; i12 < 6; i12++){
														spin.angles[12] = i12;
														for(int i13 = 0 ; i13 < 6; i13++){
															spin.angles[13] = i13;
															for(int i14 = 0 ; i14 < 6; i14++){
																spin.angles[14] = i14;
																for(int i15 = 0 ; i15 < 6; i15++){
																	spin.angles[15] = i15;
																	if(metastable(&spin, H, HB, 0)) {
																		double E = E_total(&spin, H, HB, 0) / (double)16;
																		fprintf(fichier, "%d\t%f\n", 16, E);
																		#pragma omp atomic 
																		n++;
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
  free(spin.angles);
  }
	printf("%f\t%d\n",L, n);
	free(H);
	free(HB);
	fclose(fichier);
}

void print_nxn_anti( double L, int NSimu, int nx, int ny){
	
	int Na = nx * ny + 1;
	double N = nx * ny;
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	double* n =(double*)malloc(sizeof(double) * NSimu * Na);
	for( int i = 0; i < NSimu * Na; i++) n[i] = 0.;

	
  #pragma omp parallel num_threads(10)
  {
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();

		#pragma omp for collapse(2) 
		for(int i = 0; i < Na; i++){
			for(int j = 0; j < NSimu; j++){
				spinners spin2;
				spinners_init_anti(&spin2, L, nx, ny, &seed, i);
				while (!metastable_anti(&spin2, H, HB, 0))
				{
					recuit_anti(&spin2, H, HB, 0.0001, 0.000001, 0.5, 5 * nx * ny, 0, &seed);
				}
				
				n[i * NSimu + j] = E_total_anti(&spin2, H, HB, 0) / N; 
	
    			free(spin2.angles);
				free(spin2.anti);
  			}
		}
	}


	double* nmean = (double*)malloc(Na * sizeof(double));
	double* nsd = (double*)malloc(Na * sizeof(double));

	for(int j = 0; j < Na ; j++){
		nmean[j] = 0.;
		for(int i = 0; i < NSimu; i++){
			nmean[j] += n[j * NSimu + i];
		}
		nmean[j] /= (double)NSimu;
	}

	for(int j = 0; j < Na ; j++){
		for(int i = 0; i < NSimu; i++){
			nsd[j] += (nmean[j] - n[j * NSimu + i]) * (nmean[j] - n[j * NSimu + i]);
		}
		nsd[j] = sqrt(nsd[j] / (double)NSimu);
	}


	for(int j = 0; j < Na ; j++){
		printf("%d\t%d\t%d\t%d\t%f\t%f\n",nx, ny, NSimu, j, nmean[j], nsd[j]);
	}
	
	free(nmean);
	free(nsd);
	free(n);
	free(H);
	free(HB);
	
}

void print_all3x3(char* add, double L){
	FILE* fichier = openfile_out(add);
	double* H = H_init_L(L);
	double* HB = H_B_init( 0 , 0);
	H_plot(H);
	int n = 0;
	bool *state = (bool*)malloc(6*6*6*6*6*6*6*6*6*sizeof(bool));

	int n0 = 1;
	int n1 = 6;
	int n2 = 6*6;
	int n3 = 6*6*6;
	int n4 = 6*6*6*6;
	int n5 = 6*6*6*6*6;
	int n6 = 6*6*6*6*6*6;
	int n7 = 6*6*6*6*6*6*6;
	int n8 = 6*6*6*6*6*6*6*6;
	
  #pragma omp parallel num_threads(1)
  {
  	spinners spin2;
	spinners_init(&spin2, L, 3, 3, 1);
	#pragma omp for collapse(6) 
	for(int i0 = 0 ; i0 < 6; i0++){
		for(int i1 = 0 ; i1 < 6; i1++){
			for(int i2 = 0 ; i2 < 6; i2++){
				for(int i3 = 0 ; i3 < 6; i3++){
					for(int i4 = 0 ; i4 < 6; i4++){
						for(int i5 = 0 ; i5 < 6; i5++){
							spin2.angles[0] = i0;
							spin2.angles[1] = i1;
							spin2.angles[2] = i2;
							spin2.angles[3] = i3;
							spin2.angles[4] = i4;
							spin2.angles[5] = i5;
							for(int i6 = 0 ; i6 < 6; i6++){
								spin2.angles[6] = i6;
								for(int i7 = 0 ; i7 < 6; i7++){
									spin2.angles[7] = i7;
									for(int i8 = 0 ; i8 < 6; i8++){
										spin2.angles[8] = i8;
										state[i0*n0+i1*n1+i2*n2+i3*n3+i4*n4+i5*n5+i6*n6+i7*n7+i8*n8] = metastable(&spin2, H, HB, 0);
									}
								}
							}
						}
					}
				}
			}
		}
	}
  free(spin2.angles);
  }
  
  spinners spin;
  spinners_init(&spin, L, 3, 3, 1);

  for(int i0 = 0 ; i0 < 6; i0++){
		for(int i1 = 0 ; i1 < 6; i1++){
			for(int i2 = 0 ; i2 < 6; i2++){
				for(int i3 = 0 ; i3 < 6; i3++){
					for(int i4 = 0 ; i4 < 6; i4++){
						for(int i5 = 0 ; i5 < 6; i5++){
							for(int i6 = 0 ; i6 < 6; i6++){
								for(int i7 = 0 ; i7 < 6; i7++){
									for(int i8 = 0 ; i8 < 6; i8++){
										if(state[i0*n0+i1*n1+i2*n2+i3*n3+i4*n4+i5*n5+i6*n6+i7*n7+i8*n8]) {
											spin.angles[0] = i0;
											spin.angles[1] = i1;
											spin.angles[2] = i2;
											spin.angles[3] = i3;
											spin.angles[4] = i4;
											spin.angles[5] = i5;
											spin.angles[6] = i6;
											spin.angles[7] = i7;
											spin.angles[8] = i8;	
											
											double E = E_total(&spin, H, HB, 0) / (double)9;
											fprintf(fichier, "%d\t%f\n", 16, E);
											n++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	free(spin.angles);
	free(state);
	printf("%d\t%d\n",9, n);
	free(H);
	free(HB);
	fclose(fichier);
}


void print_all3x3_anti( double L, int NSimu){
	
	double* H = H_init_L(L);
	double* HB = H_B_init( 0 , 0);
	int* n =(int*)malloc(sizeof(int) * NSimu * 10);
	double* E =(double*)malloc(sizeof(double) * NSimu * 10);
	for( int i = 0; i < NSimu * 10; i++) n[i] = 0.;

	
  #pragma omp parallel num_threads(10)
  {
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();

		#pragma omp for collapse(2) 
		for(int i = 0; i < 10; i++){
			for(int j = 0; j < NSimu; j++){
				n[i * NSimu + j] = 0;
				E[i * NSimu + j] = 0;
				spinners spin2;
				spinners_init_anti(&spin2, L, 3, 3, &seed, i);
				for(int i0 = 0 ; i0 < 6; i0++){
					for(int i1 = 0 ; i1 < 6; i1++){
						for(int i2 = 0 ; i2 < 6; i2++){
							for(int i3 = 0 ; i3 < 6; i3++){
								for(int i4 = 0 ; i4 < 6; i4++){
											for(int i5 = 0 ; i5 < 6; i5++){
										spin2.angles[0] = i0;
										spin2.angles[1] = i1;
										spin2.angles[2] = i2;
										spin2.angles[3] = i3;
										spin2.angles[4] = i4;
										spin2.angles[5] = i5;
										for(int i6 = 0 ; i6 < 6; i6++){
											spin2.angles[6] = i6;
											for(int i7 = 0 ; i7 < 6; i7++){
												spin2.angles[7] = i7;
												for(int i8 = 0 ; i8 < 6; i8++){
													spin2.angles[8] = i8;
													if(metastable_anti(&spin2, H, HB, 0)){
														#pragma omp atomic 
														n[i * NSimu + j]++;
														E[i * NSimu + j] += E_total_anti(&spin2, H, HB, 0) / 9.;
														
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
    			free(spin2.angles);
				free(spin2.anti);

				E[i * NSimu + j] /= n[i * NSimu + j];
  			}
		}
	}

	double nmean[10] = {0.};
	double nsd[10] = {0.};

	double emean[10] = {0.};
	double esd[10] = {0.};

	for(int j = 0; j < 10 ; j++){
		for(int i = 0; i < NSimu; i++){
			nmean[j] += n[j * NSimu + i];
			emean[j] += E[j * NSimu + i];
		}
		nmean[j] /= (double)NSimu;
		emean[j] /= (double)NSimu;
	}

	for(int j = 0; j < 10 ; j++){
		for(int i = 0; i < NSimu; i++){
			nsd[j] += (nmean[j] - n[j * NSimu + i]) * (nmean[j] - n[j * NSimu + i]);
			esd[j] += (emean[j] - E[j * NSimu + i]) * (emean[j] - E[j * NSimu + i]);
		}
		nsd[j] = sqrt(nsd[j] / (double)NSimu);
		esd[j] = sqrt(esd[j] / (double)NSimu);
	}


	for(int j = 0; j < 10 ; j++){
		printf("%d\t%d\t%f\t%f\t%f\t%f\n",j, NSimu, emean[j], esd[j], nmean[j], nsd[j]);
	}
	
	free(H);
	free(HB);
	free(n);
	free(E);
	
}

void plot_all3x3_of_LBX(char* add, double L, double bx, double by ,int p ){
	
	double* H = H_init_L(L);
	double* HB = H_B_init_L( bx , by, L);
	int n = 0;
	bool *state = (bool*)malloc(6*6*6*6*6*6*6*6*6*sizeof(bool));

	int n0 = 1;
	int n1 = 6;
	int n2 = 6*6;
	int n3 = 6*6*6;
	int n4 = 6*6*6*6;
	int n5 = 6*6*6*6*6;
	int n6 = 6*6*6*6*6*6;
	int n7 = 6*6*6*6*6*6*6;
	int n8 = 6*6*6*6*6*6*6*6;
	
	double M = 0;
	std::vector<double> m(0);

  #pragma omp parallel num_threads(p)
  {
  	spinners spin2;
	spinners_init(&spin2, L, 3, 3, 1);
	#pragma omp for collapse(6) 
	for(int i0 = 0 ; i0 < 6; i0++){
		for(int i1 = 0 ; i1 < 6; i1++){
			for(int i2 = 0 ; i2 < 6; i2++){
				for(int i3 = 0 ; i3 < 6; i3++){
					for(int i4 = 0 ; i4 < 6; i4++){
						for(int i5 = 0 ; i5 < 6; i5++){
							spin2.angles[0] = i0;
							spin2.angles[1] = i1;
							spin2.angles[2] = i2;
							spin2.angles[3] = i3;
							spin2.angles[4] = i4;
							spin2.angles[5] = i5;
							for(int i6 = 0 ; i6 < 6; i6++){
								spin2.angles[6] = i6;
								for(int i7 = 0 ; i7 < 6; i7++){
									spin2.angles[7] = i7;
									for(int i8 = 0 ; i8 < 6; i8++){
										spin2.angles[8] = i8;
										state[i0*n0+i1*n1+i2*n2+i3*n3+i4*n4+i5*n5+i6*n6+i7*n7+i8*n8] = metastable(&spin2, H, HB, 0);

										if(metastable(&spin2, H, HB, 0)){
											double t = magnetisation(&spin2);
											M += t;
											
											m.push_back(t);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
  free(spin2.angles);
  }
  
	double SD = 0;
	M = M / (double)m.size();

	for(int i = 0; i < m.size(); i++){
		SD += (M - m.at(i)) * (M - m.at(i));
	}

	free(state);
	printf("%f\t%f\t%d\t%f\t%f\n", bx, by,  m.size(), M, sqrt(SD/(double)m.size()));
	free(H);
	free(HB);
}


void print_all2x2(char* add, double L){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;
	
 
  	spinners spin;
	spinners_init(&spin, L, 2, 2, 1);
	for(int i0 = 0 ; i0 < 6; i0++){
		spin.angles[0] = i0;
		for(int i1 = 0 ; i1 < 6; i1++){
			spin.angles[1] = i1;
			for(int i2 = 0 ; i2 < 6; i2++){
				spin.angles[2] = i2;
				for(int i3 = 0 ; i3 < 6; i3++){
					spin.angles[3] = i3;
					if(metastable(&spin, H, HB, 0)) {
						
						double E = E_total(&spin, H, HB, 0) / (double)4;
						fprintf(fichier, "%d\t%f\n", 16, E);
						n++;
					}
				}
			}
		}
	}
	free(spin.angles);
  
	printf("%d\t%d\n",4, n);
	free(H);
	free(HB);
	fclose(fichier);
}

void print_all2x1(char* add, double L){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;
	
 
  	spinners spin;
	spinners_init(&spin, L, 2, 1, 1);
	for(int i0 = 0 ; i0 < 6; i0++){
		spin.angles[0] = i0;
		for(int i1 = 0 ; i1 < 6; i1++){
			spin.angles[1] = i1;
			if(metastable(&spin, H, HB, 0)) {
				print_spinners(&spin, fichier);
				n++;
			}
		}
	}
	free(spin.angles);
  
	printf("%f\t%d\n",L, n);
	free(H);
	free(HB);
	fclose(fichier);
}

void print_all2x1_anti(char* add, double L){
	FILE* fichier = openfile_out(add);
	double* H = H_init(L);
	double* HB = H_B_init( 0 , 0);
	int n = 0;
	
 
  	spinners spin;
	unsigned int seed = 1;
	spinners_init_anti(&spin, L, 2, 1, &seed, 1);
	spin.anti[0] = 1;
	spin.anti[1] = -1;
	for(int i0 = 0 ; i0 < 6; i0++){
		spin.angles[0] = i0;
		for(int i1 = 0 ; i1 < 6; i1++){
			spin.angles[1] = i1;
			if(metastable_anti(&spin, H, HB, 0)) {
				print_spinners(&spin, fichier);
				n++;
			}
		}
	}
	free(spin.angles);
	free(spin.anti);
  
	printf("%f\t%d\n",L, n);
	free(H);
	free(HB);
	fclose(fichier);
}

void print_Emin( spinners_t* spin,  char* add, int Niters, int p)
{
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}

	FILE* fichier = openfile_out(add);
	double* H = H_init_L(spin->L);
	double* HB = H_B_init_L(0, 0, spin->L);
	int N = spin->nx * spin->ny;
	double Emin = DBL_MAX;
	int* spinmin = (int*)malloc(N * sizeof(int));

	spin->angles = (int*)realloc(spin->angles, N * Niters * sizeof(int));
	
	double* spinE = (double*)malloc(sizeof(double) * Niters);
	bool* spinmeta = (bool*)malloc(sizeof(bool) * Niters);

	if (spinmin == NULL || spin->angles == NULL || spinE == NULL || spinmeta == NULL) {
		fprintf(stderr, "print_Emin : Allocation de memoire echouee.\n");
		exit(EXIT_FAILURE);
	}

	#pragma omp parallel num_threads(p)
	{
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();

		#pragma omp for
		for (int i = 0; i < Niters; i++) {
			for (int j = 0; j < N; j++) { spin->angles[j + N * i] = rand() % 6; }
			recuit(spin, H, HB, 1, 0.001, 0.95, 100 * N, i * N, &seed);
			recuit(spin, H, HB, 0.00000000001, 0.000000000001, 0.95, 100 * N, i * N, &seed);
			spinE[i] = E_total(spin, H, HB, N * i);
			spinmeta[i] = metastable(spin, H, HB, N * i);
		}
	}

	int iEmin = INT_MAX;
	bool nfound = true;
	for (int i = 0; i < Niters; i++) {
		if (spinE[i] < Emin && spinmeta[i]) { 
			Emin = spinE[i];
			nfound = false;
			iEmin = i;
		}
	}

	if( !nfound && iEmin != INT_MAX){
		for (int l = 0; l < N; l++) { spinmin[l] = spin->angles[l + iEmin * N] ; }
		spin->angles = (int*)realloc(spin->angles, N * sizeof(int));
		for (int l = 0; l < N; l++) { spin->angles[l] = spinmin[l] ; }
		printf("nx %d ny %d EF = %f metastable : %d\n",spin->nx, spin->ny, Emin, metastable(spin, H, HB, 0));
		print_spinners(spin, fichier);
	}
	
	free(H);
	free(HB);
	free(spinmin);
	free(spinE);
	free(spinmeta);
	fclose(fichier);
}

void print_Emax( spinners_t* spin,  char* add, int Niters){
	if(spin->Ngrid != 1){
		printf("read_spinner : invalide Ngrid %d", spin->Ngrid );
	}

	FILE* fichier = openfile_out(add);
	double* H = H_init(spin->L);
	double* HB = H_B_init(0, 0);
	int N = spin->nx * spin->ny;
	double Emax = -DBL_MAX;
	int* spinmin = (int*)malloc(N * sizeof(int));

	unsigned int seed = (unsigned int)time(NULL);

	for (int i = 0; i < Niters; i++) {
		for (int j = 0; j < N; j++) { spin->angles[j] = rand() % 6; }
		recuit(spin, H, HB, 0.0001, 0.00001, 0.5, 5 * N, 0, &seed);
		double E = E_total(spin, H, HB, 0);
		if (E > Emax && metastable(spin, H, HB, 0)) { 
			Emax = E;
			for (int l = 0; l < N; l++) { spinmin[l] = spin->angles[l]; }
		}
	}
	for (int l = 0; l < N; l++) { spin->angles[l] = spinmin[l]; }
	printf("nx = %d ny = %d, EF = %f metastable : %d\n",spin->nx, spin->ny, Emax, metastable(spin, H, HB, 0));
	print_spinners(spin, fichier);
	free(H);
	free(HB);
	free(spinmin);
	fclose(fichier);
}


void recuitN(spinners_t* spin, double* H, double* HB, double T0, double TF, double lambda, int Niter, int Nsimu, int p){

	int N = spin->nx * spin->ny;
	int* angles0 = (int*)malloc(N * sizeof(int));
	if (angles0  == NULL) {
        fprintf(stderr, "recuitN, angles0 failled to malloc\n");
        exit(EXIT_FAILURE);
	}

	memcpy(angles0, spin->angles, N * sizeof(int));
	spin->angles = (int*)realloc(spin->angles, Nsimu * N * sizeof(int));
	if (spin->angles == NULL) {
    	fprintf(stderr, "recuitN, échec de l'allocation mémoire.\n");
    	exit(EXIT_FAILURE);	
	}

	spin->Ngrid = Nsimu;
	for(int i = 0; i < spin->Ngrid; i++){
		for(int j = 0; j < N; j++){
			spin->angles[i * N + j] = angles0[j];
		}
	}
	free(angles0);

	#pragma omp parallel num_threads(p)
	{	
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();
		#pragma omp for
		for(int i = 0; i < spin->Ngrid ; i++){
			recuit(spin, H, HB, T0, TF, lambda, Niter, N * i, &seed);
		}
	}

	remove_equale_allmeta(spin, H, HB);
}

void recuitN_RSBk(spinners_t* spin, double* H, double* HB, int Nrand, double T0, double TF, double T02, double lambda, int Niter, int Nsimu, int NT02){

	int N = spin->nx * spin->ny;
	int* angles0 = (int*)malloc(N * sizeof(int) * Nrand);
	if (angles0  == NULL) {
        fprintf(stderr, "recuitN, angles0 failled to malloc\n");
        exit(EXIT_FAILURE);
	}

	spin->angles = (int*)realloc(spin->angles, Nsimu * N * sizeof(int) * Nrand * NT02);
	if (spin->angles == NULL) {
    	fprintf(stderr, "recuitN, échec de l'allocation mémoire.\n");
    	exit(EXIT_FAILURE);	
	}

	spin->Ngrid = Nsimu * Nrand * NT02;
	#pragma omp parallel
	{	
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();
		#pragma omp for
		for(int i = 0; i < Nrand; i++){
			spinners_t spin0;
			spinners_init(&spin0, spin->L, spin->nx, spin->ny, 1);
			for(int j = 0; j < N; j++) spin0.angles[j] = rand_r(&seed) % 6;
			recuit(&spin0, H, HB, 1, 0.001, 0.95, 100 * N, 0, &seed);
			recuit(&spin0, H, HB, 0.00000000001, 0.000000000001, 0.95, 100 * N, 0, &seed);
			for(int j = 0; j < N; j++) angles0[i * Nrand + j] = spin0.angles[j];
			free(spin0.angles);
		}
	}

	#pragma omp parallel for collapse(2)
	for(int i = 0; i < Nrand; i++){
		for(int j = 0; j < Nsimu * NT02; j++){
			for(int l = 0; l < N; l++){
				spin->angles[(i * Nsimu * NT02 + j) * N + l] = angles0[i * Nrand + l];
			}
		}
	}

	#pragma omp parallel
	{	
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();
		#pragma omp for
		for(int i = 0; i < Nsimu * Nrand ; i++){
			double T = T02 + (T0 - T02) * (double)rand_r(&seed) / (double)RAND_MAX;
			recuit(spin, H, HB, T, TF, lambda, Niter, N * i * NT02, &seed);
			for(int j = 0; j < NT02; j++){
				for(int l = 0; l < N; l++){
					spin->angles[N * i * NT02 + j * N + l] = spin->angles[ N * i * NT02 + l];
				}
			}
		}
	}

	#pragma omp parallel
	{	
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();
		#pragma omp for
		for(int i = 0; i < spin->Ngrid ; i++){
			recuit(spin, H, HB, T02, TF, lambda, Niter, N * i, &seed);
		}
	}

	free(angles0);

	remove_equale_allmeta(spin, H, HB);
}

void recuitN_RSBk_anti(spinners_t* spin, double* H, double* HB, int Nrand, double T0, double TF, double T02, double lambda, int Niter, int Nsimu, int NT02, unsigned int* seed, double prob_anti){

	int N = spin->nx * spin->ny;
	int* angles0 = (int*)malloc(N * sizeof(int) * Nrand);
	if (angles0  == NULL) {
        fprintf(stderr, "recuitN, angles0 failled to malloc\n");
        exit(EXIT_FAILURE);
	}

	int* angles0_anti = (int*)malloc(N * sizeof(int) * Nrand);
	if (angles0_anti  == NULL) {
        fprintf(stderr, "recuitN, angles0_anti failled to malloc\n");
        exit(EXIT_FAILURE);
	}

	spin->angles = (int*)realloc(spin->angles, Nsimu * N * sizeof(int) * Nrand * NT02);
	if (spin->angles == NULL) {
    	fprintf(stderr, "recuitN, échec de l'allocation mémoire.\n");
    	exit(EXIT_FAILURE);	
	}

	spin->anti = (int*)realloc(spin->anti, Nsimu * N * sizeof(int) * Nrand * NT02);
	if (spin->angles == NULL) {
    	fprintf(stderr, "recuitN, échec de l'allocation mémoire.\n");
    	exit(EXIT_FAILURE);	
	}

	spin->Ngrid = Nsimu * Nrand * NT02;
	for(int i = 0; i < Nrand; i++){
			spinners_t spin0;
			spinners_init_anti_prob(&spin0, spin->L, spin->nx, spin->ny, 1, seed, prob_anti);
			for(int j = 0; j < N; j++) spin0.angles[j] = rand_r(seed) % 6;
			recuit_anti(&spin0, H, HB, 1, 0.001, 0.95, 100 * N, 0, seed);
			recuit_anti(&spin0, H, HB, 0.00000000001, 0.000000000001, 0.95, 100 * N, 0, seed);
			for(int j = 0; j < N; j++) {
				angles0[i * Nrand + j] = spin0.angles[j];
				angles0_anti[i * Nrand + j] = spin0.anti[j];
			}
			free(spin0.angles);
		}


	for(int i = 0; i < Nrand; i++){
		for(int j = 0; j < Nsimu * NT02; j++){
			for(int l = 0; l < N; l++){
				spin->angles[(i * Nsimu * NT02 + j) * N + l] = angles0[i * Nrand + l];
				spin->anti[(i * Nsimu * NT02 + j) * N + l] = angles0_anti[i * Nrand + l];
			}
		}
	}

	for(int i = 0; i < Nsimu * Nrand ; i++){
		double T = T02 + (T0 - T02) * (double)rand_r(seed) / (double)RAND_MAX;
		recuit_anti(spin, H, HB, T, TF, lambda, Niter, N * i * NT02, seed);
		for(int j = 0; j < NT02; j++){
			for(int l = 0; l < N; l++){
				spin->angles[N * i * NT02 + j * N + l] = spin->angles[ N * i * NT02 + l];
				spin->anti[N * i * NT02 + j * N + l] = spin->anti[ N * i * NT02 + l];
			}
		}
	}

	for(int i = 0; i < spin->Ngrid ; i++){
			recuit_anti(spin, H, HB, T02, TF, lambda, Niter, N * i, seed);
		}

	free(angles0);
	free(angles0_anti);

	remove_equale_allmeta_anti(spin, H, HB);
}


void UM_of_L(int nx, int ny, double L, double T0, double T02, int Nrand, int Niter, int Nmean, int NT02){

	if (T0 <= 0.00001) {
    	fprintf(stderr, "T0 is to low\n");
    	exit(EXIT_FAILURE);	
	}

	double* dEG = (double*)malloc(sizeof(double) * Nmean);
	double* dEL = (double*)malloc(sizeof(double) * Nmean);
	double* dH = (double*)malloc(sizeof(double) * Nmean);
	double* dIH = (double*)malloc(sizeof(double) * Nmean);

	for(int i = 0; i < Nmean; i++){
		spinners_t spin;
		spinners_init(&spin, L, nx, ny, 1);

		double* H = H_init_L25(L);
		double* HB = H_B_init_L(0, 0, L);

		double degree[4] = {0};
		recuitN_RSBk(&spin, H, HB, Nrand, T0, 0.00001, T02, 0.1, nx * nx, Niter, NT02);
		plot_degree_UM_and_stat(&spin, degree);

		dEG[i] = degree[0];
		dEL[i] = degree[1];
		dH[i] = degree[2];
		dIH[i] = degree[3];
	
		free(H);
		free(HB);
		free(spin.angles);
	}

	printf("%f\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", L, nx, mean(dEG, Nmean), SD(dEG, Nmean), mean(dEL, Nmean), SD(dEL, Nmean)
										, mean(dH, Nmean), SD(dH, Nmean), mean(dIH, Nmean), SD(dIH, Nmean));
	fflush(stdout);

	free(dEG);
	free(dEL);
	free(dH);
	free(dIH);
}

void UM_of_L_anti(int nx, int ny, double L, double T0, double T02, int Nrand, int Niter, int NT02, double prob_anti, int Nquechedrandom){

	if (T0 <= 0.00001) {
    	fprintf(stderr, "T0 is to low\n");
    	exit(EXIT_FAILURE);	
	}

	double* dEL = (double*)malloc(sizeof(double)  * Nquechedrandom);

	#pragma omp parallel
	{

		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();

		#pragma omp for
		for(int h = 0; h < Nquechedrandom; h++){
			
				spinners_t spin;
				spinners_init_anti_prob(&spin, L, nx, ny, 1, &seed, prob_anti);

				double* H = H_init_L25(L);
				double* HB = H_B_init_L(0, 0, L);

				double degree;
				recuitN_RSBk_anti(&spin, H, HB, Nrand, T0, 0.00001, T02, 0.1, nx * nx, Niter, NT02, &seed, prob_anti);
				plot_degree_UM_and_stat_anti(&spin, &degree);

				dEL[h] = degree;
	
				free(H);
				free(HB);
				free(spin.angles);
				free(spin.anti);
			
		}
	}



	printf("%f\t%f\t%d\t%d\t%f\t%f\n", L, prob_anti, nx, ny, mean(dEL, Nquechedrandom), sqrt(SD(dEL, Nquechedrandom)));
	fflush(stdout);

	free(dEL);
}


void insert_ones(int index, int remaining, int configurations[512][9], int *count, int a[9]) {
	if (remaining == 0) {
		
		for (int i = 0; i < 9; i++) {
			configurations[*count][i] = a[i];
		}
		(*count)++; 
		return;
	}

	
	for (int i = index; i < 9; i++) {
		a[i] = -1; 
		insert_ones(i + 1, remaining - 1, configurations, count, a);
		a[i] = 1; 
	}
}

void generate_configurations(int n, int configurations[512][9], int *count) {
    
    if (n < 0 || n > 9) {
        printf("Le nombre de 1 doit être compris entre 0 et 9.\n");
        return;
    }

    int a[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}; 

   
    insert_ones(0, n, configurations, count, a);
}

void print_configurations(int configurations[512][9], int count) {
    for (int i = 0; i < count; i++) {
        printf("[");
        for (int j = 0; j < 9; j++) {
            printf("%d", configurations[i][j]);
            if (j < 8) printf(", ");
        }
        printf("]\n");
    }
	printf("count %d\n", count);
}

void anti9_all(double L){
	
	int NSimu = 126;
	double* H = H_init_L(L);
	double* HB = H_B_init( 0 , 0);
	int* n =(int*)malloc(sizeof(int) * NSimu * 10);
	double* E =(double*)malloc(sizeof(double) * NSimu * 10);
	for( int i = 0; i < NSimu * 10; i++) n[i] = 0.;

	int NSimumax[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	
  #pragma omp parallel num_threads(10)
  {
		unsigned int seed = (unsigned int)time(NULL) + omp_get_thread_num();

		#pragma omp for collapse(1) 
		for(int i = 0; i < 10; i++){

			int configurations[512][9];

			generate_configurations(i, configurations, &NSimumax[i]);
			
			for(int j = 0; j < NSimumax[i]; j++){
				printf("%d\t%d\t%d\n", NSimumax[i], i, j);
				n[i * NSimu + j] = 0;
				E[i * NSimu + j] = 0;
				spinners spin2;
				spinners_init_anti(&spin2, L, 9, 1, &seed, 0);
				
				spin2.anti[0] = configurations[j][0];
				spin2.anti[1] = configurations[j][1];
				spin2.anti[2] = configurations[j][2];
				spin2.anti[3] = configurations[j][3];
				spin2.anti[4] = configurations[j][4];
				spin2.anti[5] = configurations[j][5];
				spin2.anti[6] = configurations[j][6];
				spin2.anti[7] = configurations[j][7];
				spin2.anti[8] = configurations[j][8];

				for(int i0 = 0 ; i0 < 6; i0++){
					for(int i1 = 0 ; i1 < 6; i1++){
						for(int i2 = 0 ; i2 < 6; i2++){
							for(int i3 = 0 ; i3 < 6; i3++){
								for(int i4 = 0 ; i4 < 6; i4++){
											for(int i5 = 0 ; i5 < 6; i5++){
										spin2.angles[0] = i0;
										spin2.angles[1] = i1;
										spin2.angles[2] = i2;
										spin2.angles[3] = i3;
										spin2.angles[4] = i4;
										spin2.angles[5] = i5;
										for(int i6 = 0 ; i6 < 6; i6++){
											spin2.angles[6] = i6;
											for(int i7 = 0 ; i7 < 6; i7++){
												spin2.angles[7] = i7;
												for(int i8 = 0 ; i8 < 6; i8++){
													spin2.angles[8] = i8;
													if(metastable_anti(&spin2, H, HB, 0)){
														#pragma omp atomic 
														n[i * NSimu + j]++;
														E[i * NSimu + j] += E_total_anti(&spin2, H, HB, 0) / 9.;
														
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
    			free(spin2.angles);
				free(spin2.anti);

				E[i * NSimu + j] /= n[i * NSimu + j];
  			}
		}
	}

	double nmean[10] = {0.};
	double nsd[10] = {0.};

	double emean[10] = {0.};
	double esd[10] = {0.};

	for(int j = 0; j < 10 ; j++){
		for(int i = 0; i < NSimumax[j]; i++){
			nmean[j] += n[j * NSimu + i];
			emean[j] += E[j * NSimu + i];
		}
		nmean[j] /= (double)NSimumax[j];
		emean[j] /= (double)NSimumax[j];
	}

	for(int j = 0; j < 10 ; j++){
		for(int i = 0; i < NSimumax[j]; i++){
			nsd[j] += (nmean[j] - n[j * NSimu + i]) * (nmean[j] - n[j * NSimu + i]);
			esd[j] += (emean[j] - E[j * NSimu + i]) * (emean[j] - E[j * NSimu + i]);
		}
		nsd[j] = sqrt(nsd[j] / (double)NSimu);
		esd[j] = sqrt(esd[j] / (double)NSimu);
	}


	for(int j = 0; j < 10 ; j++){
		printf("%d\t%d\t%f\t%f\t%f\t%f\n",j, NSimumax[j], emean[j], esd[j], nmean[j], nsd[j]);
	}
	
	free(H);
	free(HB);
	free(n);
	free(E);

}



