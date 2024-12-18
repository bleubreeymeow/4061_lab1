#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define INITIAL_STR "initial_coords"
#define PERTURBE_STR "perturbed_coords"
#define NEIGH_LIST_STR "neighb_list"
#define FINAL_STR "final"
#define SIGMA 3.40 //sigma for Ar (unit of angstroms) 
#define EPSILON 0.0104233206
#define UX 5
#define UY 5
#define UZ 5
#define A_CONSTANT 5.31


double distance_cutoff = 2.5 * SIGMA;  
int atom_num = UX * UY * UZ;
int fcc_atom_num = UX * UY * UZ * 4;
double volume;
double b[3][3]; //reciprocal unit cell vectors
double a[3][3]; //unit cell vectors

void FILE_WRITING(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "LAB4_1_%s_SDLJ.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%lf \t %lf \t %lf \n",arr[i][0], arr[i][1] , arr[i][2]);
    }
    printf("%s has been printed!\n",lattice_structure);
    fclose(filepointer);
}

void LJ_FILE_WRITING(double *arr, int num){    
    FILE *filepointer = fopen("LAB4_1_SDLJ_energy_vs_steps.txt", "w");
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%lf \n",arr[i]);
    }
    printf("energy vs steps has been printed!\n");
    fclose(filepointer);
}


/*b1 b2 b3 VECTOR GENERATION=============================================================================================*/
void fn_a_vect(){
    //generating a1 a2 a3 vectors
    a[0][0] = UX * A_CONSTANT;
    a[1][1] = UY * A_CONSTANT;
    a[2][2] = UZ * A_CONSTANT;
    for(int i = 0 ; i < 3 ; i++){
        for(int j = 0 ; j < 3 ; j++){
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

void fn_reciprocal(){
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

/*ATOM COORDS GENERATION==============================================================================================*/

//shifting coords, for translation
void fn_shift_coords(double **coords, int count, double shift_factor){
    for (int i = 0; i < count; i++) {
        for (int j = 0; j < 3; j++) {
            coords[i + count][j] = coords[i][j] + A_CONSTANT * shift_factor;
        }
    }
    return;
}

void fn_simplecubic(double **coords){
    for (int i = 0; i < atom_num ; i++) {
        coords[i][0] = (i % UX) * A_CONSTANT;
        coords[i][1] = (i / UX) % UY * A_CONSTANT;
        coords[i][2] = (i / UX / UY) % UZ * A_CONSTANT;
    }
    return;
}

void fn_fcc(double **coords){
    fn_simplecubic(coords); //construct a simple cubic lattice
    for (int i = atom_num; i < atom_num * 4; i += 3) {
        int index = (int)((i - atom_num) / 3); //retrieve the atom index in the simple cubic lattice
        //2nd basis atom at (a/2 , a/2 , 0) relative to 1st basis atom
        coords[i][0] = coords[index][0] + A_CONSTANT * 0.5;
        coords[i][1] = coords[index][1] + A_CONSTANT * 0.5;
        coords[i][2] = coords[index][2];
        //3rd basis atom at (a/2 , 0 , a/2) relative to 1st basis atom       
        coords[i + 1][0] = coords[index][0] + A_CONSTANT * 0.5;
        coords[i + 1][1] = coords[index][1];
        coords[i + 1][2] = coords[index][2] + A_CONSTANT * 0.5;
        //4th basis atom at (0 , a/2 , a/2) relative to 1st basis atom
        coords[i + 2][0] = coords[index][0];
        coords[i + 2][1] = coords[index][1] + A_CONSTANT * 0.5;
        coords[i + 2][2] = coords[index][2] + A_CONSTANT * 0.5;
    }
    return;
}


/*INITIALISATION=====================================================================================================*/
void FN_initialise(double** atom_coords){

    fn_a_vect(); //initialise a1 a2 a3 cell vectors
    fn_reciprocal(); //calculate the reciprocal vectors and the volume
    fn_fcc(atom_coords); //generate fcc coordinates

    return;
}


/*PERTURBATING ATOMS================================================================================================*/
void FN_perturbate(double** atom_coords){
    double min = 0.0;
    double max = 0.5;
    for(int i = 0 ; i < fcc_atom_num ; i++){
        for(int j = 0 ; j < 3 ; j++){
            double random = (rand() * (max - min) / RAND_MAX ) + min;
            atom_coords[i][j] += random;
        }
    }

    return;
}

/*NEIGHBOUR LIST GENERATION=======================================================================================*/

double fn_pbc(double* t){
    double n[3];
    double difference_plus = 0 , difference_minus = 0;
    for(int i = 0 ; i < 3 ; i++){
        //generate fractional numbers
        n[i] = fmod((t[0] * b[i][0] + t[1] * b[i][1] + t[2] * b[i][2]) , 1);
        //apply periodic boundary conditions [-0.5 , 0.5)
        if(n[i] > 0.5){
            difference_plus = n[i] - 0.5;
            n[i] = - 0.5 + difference_plus;
            difference_plus = 0;
        }
        if(n[i] <= -0.5){
            difference_minus = -n[i] - 0.5;
            n[i] = 0.5 - difference_minus;
            difference_minus = 0;
        }
    }
    for(int i = 0 ; i < 3 ; i++){
        t[i] = n[0] * a[0][i] + n[1] * a[1][i] + n[2] * a[2][i]; //fractional coordinate
    }
    double distance = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
    return distance;
}

int fn_neighbour_list(double **neigh_list, double **atom_coords, int size){
    double distance = 0;
    int nearest_neighbour_num = 0;
    double t[3];

    for(int i = 0 ; i < size ; i++){
        for(int j = i + 1 ; j < size ; j++){ //j = i + 1 to avoid double counting (1 2 , 2 1), i + 1 for skipping (0 0) , (1 1)

            t[0] = atom_coords[j][0] - atom_coords[i][0]; //dx
            t[1] = atom_coords[j][1] - atom_coords[i][1]; //dy
            t[2] = atom_coords[j][2] - atom_coords[i][2]; //dz
            distance = fn_pbc(t);

            if(distance < distance_cutoff){ //evaulate whether atom i and atom j are nearest neighbours

                neigh_list[nearest_neighbour_num][0] = i;
                neigh_list[nearest_neighbour_num][1] = j;
                neigh_list[nearest_neighbour_num][2] = distance;
                nearest_neighbour_num = nearest_neighbour_num + 1;
            }
        }
    }

    return nearest_neighbour_num;
}


/*LJ POTENTIAL CALCULATION==========================================================================================*/

double fn_LJ_potential(double** neighbour_list, int size){
    double potential = 0;
    for(int i = 0 ; i < size ; i++){
        double R_3 = (SIGMA / neighbour_list[i][2]) * (SIGMA / neighbour_list[i][2]) * (SIGMA / neighbour_list[i][2]);
        double R_6 = R_3 * R_3;
        potential += 4* EPSILON *((R_6 * R_6) - R_6);
    }
    return potential;
}


void FN_SD(double** atom_coords){
    int i_max = 100;
    int i = 0;

    //create neighbour list
    double row_num = fcc_atom_num * (fcc_atom_num - 1);
    double **neighbour_list = (double **)malloc(row_num * sizeof(double *));
    for (int j = 0; j < row_num; j++){
        neighbour_list[j] = (double *)malloc(3 * sizeof(double));
    }

    //create array to hold energy for each SD step
    double *LJ_energy;
    LJ_energy = (double*)malloc(i_max * sizeof(double));

    int neighbour_num = fn_neighbour_list(neighbour_list, atom_coords, fcc_atom_num);
    FILE_WRITING(NEIGH_LIST_STR,neighbour_list,neighbour_num);

    LJ_energy[0] = fn_LJ_potential(neighbour_list, neighbour_num);
    printf("LJ energy per atom (eV) = %lf \n",LJ_energy[0]/fcc_atom_num);

    for(i = 0 ; i < i_max ; i++){
        


    }



    LJ_FILE_WRITING(LJ_energy,i);
    free(LJ_energy);
    return;
}




/*MAIN==============================================================================================================*/
int main(){
    srand(time(NULL));
    double** fcc_atom_coords = (double **)malloc(fcc_atom_num * sizeof(double *));
    for (int i = 0; i < fcc_atom_num; i++){
        fcc_atom_coords[i] = (double *)malloc(3 * sizeof(double));
    }

    FN_initialise(fcc_atom_coords);
    //initial unpreturbed atom coords
    FILE_WRITING(INITIAL_STR , fcc_atom_coords,fcc_atom_num);

    //perturbed atom coords
    FN_perturbate(fcc_atom_coords);
    FILE_WRITING(PERTURBE_STR , fcc_atom_coords,fcc_atom_num);


    //perform steepest descent
    FN_SD(fcc_atom_coords);

    free(fcc_atom_coords);
    
    return 0;
}