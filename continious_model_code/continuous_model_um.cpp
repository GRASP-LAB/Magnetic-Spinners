#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "omp.h"
#include "continuous_model_um.h"
#include "continuous_model_function.h"


t_spinner_grid* generate_state( int nx, int ny, double L, unsigned int* seed) {
    
    t_spinner_grid grid0;
    init_grid(&grid0, nx, ny, L);
    init_rand(&grid0, seed);

    annealing(&grid0, 1, 0.001, 0.5, 1000, seed);

    t_spinner_grid* gridfirst = (t_spinner_grid*)malloc(sizeof(t_spinner_grid) * 10);

    for(int i = 0; i < 10; i++){
        init_grid(&gridfirst[i], nx, ny, L);
        for(int j = 0; j < nx * ny; j++) {
            gridfirst[i].spin[j].theta = grid0.spin[j].theta; 
        }
        double ti = (0.5 - 0.04) * rand_r(seed) / (double)RAND_MAX + 0.04;
        annealing(&gridfirst[i], ti, 0.001, 0.95, 1000, seed);
    }
    
    t_spinner_grid* grid = (t_spinner_grid*)malloc(sizeof(t_spinner_grid) * 128 * 10);
    for(int i = 0; i < 10; i++){
        for(int l = 0; l < 128; l++){
            init_grid(&grid[i * 128 + l], nx, ny, L);
            for(int j = 0; j < nx * ny; j++) {
                grid[i * 128 + l].spin[j].theta = gridfirst[i].spin[j].theta; 
            }
            annealing(&grid[i * 128 + l], 0.04, 0.001, 0.95, 1000, seed);
        }
    }

    free(grid0.spin);
    for(int i = 0; i < 10; i++) free(gridfirst[i].spin);
    free(gridfirst);

    return grid;
}


double dist_EL(t_spinner_grid* grid, int i, int j){
	double d = 0.;
	for(int k = 0; k < grid[i].Nspin; k++){
        double E = V_local(&grid[i], k) - V_local(&grid[j], k);
        d += E *E ;
    }
	return std::sqrt(d / (double)grid[i].Nspin);
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

double degree_UM(t_spinner_grid* grid, int Ngrid, double (*dist)(t_spinner_grid* grid, int i, int j), 
				double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist ) ){
	double sum_diff = 0;
	double sum_metric = 0;

	matrice_t matrice;
	matrice_t ultra;
	matrice.N = Ngrid;
	ultra.N = Ngrid;
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

	for(int i = 0; i < Ngrid ; i ++){
		for(int j = i; j < Ngrid ; j ++){
			double d = dist(grid, i, j);
			matrice.line[i].col[j] = d;
			matrice.line[j].col[i] = d;
		}
	}

	matrice_ultra(&matrice, &ultra, similarity);
	
	for(int i = 0; i < Ngrid ; i ++){
		for(int j = 0; j < i ; j ++){
			sum_diff += matrice.line[i].col[j] - ultra.line[i].col[j];
			sum_metric += matrice.line[i].col[j];
		}
	}

	for(int i = 0; i < matrice.N; i++){free(matrice.line[i].col);}
	free(matrice.line);
	for(int i = 0; i < ultra.N; i++){free(ultra.line[i].col);}
	free(ultra.line);
	return abs(sum_diff)/ sum_metric;
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

FILE* openfile_out(char* add) {

	 FILE* fichier = fopen(add, "w");

    if (fichier == NULL) {
        fprintf(stderr, "openfile_out : Impossible d'ouvrir le fichier %s pour l'�criture.\n", add);
        exit(EXIT_FAILURE);
    }
	return fichier;
}

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


void print_dist(t_spinner_grid* grid, int Ngrid, char* add, double (*dist)(t_spinner_grid* grid, int i, int j), 
		double (similarity)(cluster_t cluster_A, cluster_t cluster_B, matrice_t* matrice_dist )){
	const int STRING_MAX = 256;

	char path[STRING_MAX];
	strcpy(path, add);
	strcat(path, "_L"); 
	snprintf(path + strlen(path), sizeof(path) - strlen(path), "%f", grid->L);
	strcat(path, "_distEL"); 
	strcat(path, ".txt");


	char path2[STRING_MAX];
	strcpy(path2, add);
	strcat(path2, "_L"); 
	snprintf(path2 + strlen(path2), sizeof(path2) - strlen(path2), "%f", grid->L);
	strcat(path2, "_distLE"); 
	strcat(path2, "_ultra.txt");

	FILE* fichier = openfile_out(path);
	FILE* fichier_ultra = openfile_out(path2);

	matrice_t matrice;
	matrice_t ultra;
	matrice.N = Ngrid;
	ultra.N = Ngrid;
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

	const int N = grid->Nspin;
	for(int i = 0; i < Ngrid ; i ++){
		for(int j = i; j < Ngrid ; j ++){
			double d = dist(grid, i, j);
			matrice.line[i].col[j] = d;
			matrice.line[j].col[i] = d;
		}
	}
	
	matrice_ultra(&matrice, &ultra, similarity);
	tri(&matrice, &ultra);
	print_matrice(&matrice, fichier);
	print_matrice(&ultra, fichier_ultra);
	fclose(fichier);
	fclose(fichier_ultra);
	for(int i = 0; i < matrice.N; i++){free(matrice.line[i].col);}
	free(matrice.line);
	for(int i = 0; i < ultra.N; i++){free(ultra.line[i].col);}
	free(ultra.line);
}