#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#include <string.h>

double distance_cutoff = 10* (5.64); //large enough cutoff
double a_constant = 5.25;
int ux = 5 , uy = 5;
int atom_num = ux * uy;

double b[2][2]; //reciprocal unit cell vectors
double a[2][2]; //unit cell vectors
double volume;

double **atom_coords = NULL;
double **neighbour_list = NULL;

void file_writing(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "ising_%s_neighb_list.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%.0lf \t %.0lf \t %lf \n",arr[i][0], arr[i][1] , arr[i][2]);
    }
    printf("neighbour list of %s has been printed!\n",lattice_structure);
    fclose(filepointer);
}

void FN_a_vect(){
    //generating a1 a2 a3 vectors
    a[0][0] = ux * a_constant;
    a[1][1] = uy * a_constant;
    for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 2 ; j++){
            if(i != j){
                a[i][j] = 0;
            }
        }
    }
    return;
}

void fn_cross_product(int idx1 , int idx2 , int idx3){ //b_(idx1 + 1) = a_(idx2 + 1) x a_(idx3 + 1)
        b[idx1][0] = a[idx2][1]*a[idx3][2] - a[idx2][2]*a[idx3][1];
        b[idx1][1] = a[idx2][2]*a[idx3][0] - a[idx2][0]*a[idx3][2];
        b[idx1][2] = a[idx2][0]*a[idx3][1] - a[idx2][1]*a[idx3][0];
    return;
}

double fn_2pi_triple_product(){
    volume = (a[0][0] * b[0][0] + a[0][1] * b[0][1] + a[0][2] * b[0][2]);
    double answer = 1/volume;
    return answer;
}

void FN_reciprocal(){
    fn_cross_product(0 , 1 , 2); // b1 = a2 x a3
    double scalars = fn_2pi_triple_product(); //a1 dot (a2 x a3)
    fn_cross_product(1 , 2 , 0); // b2 = a3 x a1
    fn_cross_product(2 , 0 , 1); // b3 = a1 x a2

    for (int i = 0 ; i < 3 ; i++){
        for(int j = 0 ; j < 3 ; j++){
            b[i][j] = b[i][j] * scalars;
        }
    }
    return;
}

void fn_simplecubic(double **atom_coords){
    int atom_type = 1;
    for (int i = 0; i < atom_num ; i++) {
        atom_coords[i][0] = (i % ux) * a_constant;
        atom_coords[i][1] = (i / ux) % uy * a_constant;

        atom_type *= -1;
        if(ux % 2 == 0 && i % (ux * uy) == 0){
            atom_type *= -1;
        }
        if(ux % 2 == 0 && i % ux == 0){ //if perodicity of y is an even number and the loop has reached to the end, then the atom in the next y 
            atom_type *= -1;
        }
        atom_coords[i][2] = atom_type; //all the 'simple cubic' lattice points are Na
  
    }
    return;
}

void FN_atom_coords(){


    //initialise 2D array for atom coords (initialise a atom_num x 2 array)
    atom_coords = (double **)malloc(atom_num* sizeof(double *));
    for (int i = 0; i < atom_num; i++){
        atom_coords[i] = (double *)malloc(3 * sizeof(double));
    }

    fn_simplecubic(atom_coords);

    return;
}

int fn_neighbour_list(double **target, double **source, int size, double distance_cutoff){
    double distance = 0;
    int nearest_neighbour_num = 0;
    double t[3];

    for(int i = 0 ; i < size ; i++){
        for(int j = 0 ; j < size ; j++){ //j = i + 1 to avoid double counting (1 2 , 2 1), i + 1 for skipping (0 0) , (1 1)
            if(i == j){continue;}

            t[0] = source[j][0] - source[i][0]; //dx
            t[1] = source[j][1] - source[i][1]; //dy
            t[2] = source[j][2] - source[i][2]; //dz
            distance = fn_pbc(t);

            if(distance <= distance_cutoff){ //evaulate whether atom i and atom j are nearest neighbours
                target[nearest_neighbour_num][0] = i;
                target[nearest_neighbour_num][1] = j;
                target[nearest_neighbour_num][2] = distance;
                nearest_neighbour_num = nearest_neighbour_num + 1;
            }
        }
    }
    return nearest_neighbour_num;
}


void FN_generate_neighbour_list(char *name, double** coords, int size){
    double row_num = size * (size - 1);
    double **neighbour_list = (double **)malloc(row_num * sizeof(double *));
    for (int i = 0; i < row_num; i++){
        neighbour_list[i] = (double *)malloc(3 * sizeof(double));
    }

    int neighbour_num = fn_neighbour_list(neighbour_list, coords, size);
    double LJ_potential = fn_LJ_potential(neighbour_list, neighbour_num);

    printf("Ar LJ potential (eV): %lf\n", (LJ_potential/(size)));

    free(neighbour_list);
    free(coords);
}

int main(){
    FN_a_vect(); //initialise a1 a2 a3 cell vectors
    FN_reciprocal(); //calculate the reciprocal vectors and the volume

    FN_atom_coords(); //generate atom coordinates


    return 0;
}