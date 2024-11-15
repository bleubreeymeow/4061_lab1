#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define FCC_STR "fcc3"
#define SIGMA 3.40 //sigma for Ar (unit of angstroms) 
#define EPSILON 0.0104233206

int ux = 5 , uy = 5 , uz = 5;
double a_constant = 5.31;
double distance_cutoff = 2.5 * SIGMA;
int atom_num;
int fcc_atom_num;
double volume;
double b[3][3]; //reciprocal unit cell vectors
double a[3][3]; //unit cell vectors

double **fcc_atom_coords = NULL;
double **fcc_neighbour_list = NULL;


void FN_a_vect();
void fn_cross_product(int idx1 , int idx2 , int idx3);
double fn_2pi_triple_product();
void FN_reciprocal();
void fn_shift_coords(double **atom_coords, int count, double shift_factor);
void fn_simplecubic(double **atom_coords);
double fn_LJ_potential(double** neighbour_list, int size);

/*b1 b2 b3 VECTOR GENERATION=============================================================================================*/
void FN_a_vect(){
    //generating a1 a2 a3 vectors
    a[0][0] = ux * a_constant;
    a[1][1] = uy * a_constant;
    a[2][2] = uz * a_constant;
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
/*ATOM COORDS GENERATION==============================================================================================*/

//shifting coords, for translation
void fn_shift_coords(double **atom_coords, int count, double shift_factor){
    for (int i = 0; i < count; i++) {
        for (int j = 0; j < 3; j++) {
            atom_coords[i + count][j] = atom_coords[i][j] + a_constant * shift_factor;
        }
    }
}

void fn_simplecubic(double **atom_coords){
    for (int i = 0; i < atom_num ; i++) {
        atom_coords[i][0] = (i % ux) * a_constant;
        atom_coords[i][1] = (i / ux) % uy * a_constant;
        atom_coords[i][2] = (i / ux / uy) % uz * a_constant;
    }
    return;
}

void fn_fcc(double **atom_coords){
    fn_simplecubic(atom_coords); //construct a simple cubic lattice
    for (int i = atom_num; i < atom_num * 4; i += 3) {
        int index = (int)((i - atom_num) / 3); //retrieve the atom index in the simple cubic lattice
        //2nd basis atom at (a/2 , a/2 , 0) relative to 1st basis atom
        atom_coords[i][0] = atom_coords[index][0] + a_constant * 0.5;
        atom_coords[i][1] = atom_coords[index][1] + a_constant * 0.5;
        atom_coords[i][2] = atom_coords[index][2];
        //3rd basis atom at (a/2 , 0 , a/2) relative to 1st basis atom       
        atom_coords[i + 1][0] = atom_coords[index][0] + a_constant * 0.5;
        atom_coords[i + 1][1] = atom_coords[index][1];
        atom_coords[i + 1][2] = atom_coords[index][2] + a_constant * 0.5;
        //4th basis atom at (0 , a/2 , a/2) relative to 1st basis atom
        atom_coords[i + 2][0] = atom_coords[index][0];
        atom_coords[i + 2][1] = atom_coords[index][1] + a_constant * 0.5;
        atom_coords[i + 2][2] = atom_coords[index][2] + a_constant * 0.5;
    }
    return;
}


void FN_atom_coords(){
    atom_num = ux * uy * uz;
    fcc_atom_num = atom_num * 4;

    //initialise 2D array for atom coords (initialise a atom_num x 3 array)
    fcc_atom_coords = (double **)malloc(fcc_atom_num* sizeof(double *));
    for (int i = 0; i < fcc_atom_num; i++){
        fcc_atom_coords[i] = (double *)malloc(3 * sizeof(double));
    }

    fn_fcc(fcc_atom_coords);

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

int fn_neighbour_list(double **target, double **source, int size){
    double distance = 0;
    int nearest_neighbour_num = 0;
    double t[3];

    for(int i = 0 ; i < size ; i++){
        for(int j = i + 1 ; j < size ; j++){ //j = i + 1 to avoid double counting (1 2 , 2 1), i + 1 for skipping (0 0) , (1 1)

            t[0] = source[j][0] - source[i][0]; //dx
            t[1] = source[j][1] - source[i][1]; //dy
            t[2] = source[j][2] - source[i][2]; //dz
            distance = fn_pbc(t);

            if(distance < distance_cutoff){ //evaulate whether atom i and atom j are nearest neighbours

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

/*MAIN==============================================================================================================*/
int main(){

    FN_a_vect(); //initialise a1 a2 a3 cell vectors
    FN_reciprocal(); //calculate the reciprocal vectors and the volume
    FN_atom_coords(); //generate sc, bcc, fcc and diamond coordinates

    //generate fcc neighbour list
    FN_generate_neighbour_list(FCC_STR, fcc_atom_coords, fcc_atom_num);

    return 0;
}

