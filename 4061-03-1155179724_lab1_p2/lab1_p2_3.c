#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SC_STR "sc2"
#define BCC_STR "bcc2"
#define FCC_STR "fcc2"
#define DIAMOND_FCC_STR "diamond2"

int ux , uy , uz;
double a_constant;
int atom_num;
int bcc_atom_num;
int fcc_atom_num;
int diamond_atom_num;
double volume;
double b[3][3]; //reciprocal unit cell vectors
double a[3][3]; //unit cell vectors

double **sc_atom_coords = NULL;
double **bcc_atom_coords = NULL;
double **fcc_atom_coords = NULL;
double **diamond_atom_coords = NULL;

double **sc_neighbour_list = NULL;
double **bcc_neighbour_list = NULL;
double **fcc_neighbour_list = NULL;
double **diamond_neighbour_list = NULL;

void read_input();
void file_writing(char *lattice_structure, double **arr, int num);

void FN_a_vect();
void fn_cross_product(int idx1 , int idx2 , int idx3);
double fn_2pi_triple_product();
void FN_reciprocal();

void FN_atom_coords();
void fn_shift_coords(double **atom_coords, int count, double shift_factor);
void fn_simplecubic(double **atom_coords);
void fn_bcc(double **atom_coords);
void fn_fcc(double **atom_coords);
void fn_diamond_fcc(double **atom_coords);

void FN_generate_neighbour_list(char *name, double** coords, int size);
double fn_pbc(double* t);
int fn_neighbour_list(double **target, double **source, int size, double distance_cutoff);
double fn_distance_cutoff(char *name);

/*MAIN==============================================================================================================*/
int main(){
    read_input(); 

    FN_a_vect(); //initialise a1 a2 a3 cell vectors
    FN_reciprocal(); //calculate the reciprocal vectors and the volume
    FN_atom_coords(); //generate sc, bcc, fcc and diamond coordinates

    //generate sc neighbour list
    FN_generate_neighbour_list(SC_STR, sc_atom_coords, atom_num);
    //generate bcc neighbour list
    FN_generate_neighbour_list(BCC_STR, bcc_atom_coords, bcc_atom_num);
    //generate fcc neighbour list
    FN_generate_neighbour_list(FCC_STR, fcc_atom_coords, fcc_atom_num);
    //generate diamond neighbour list
    FN_generate_neighbour_list(DIAMOND_FCC_STR, diamond_atom_coords, diamond_atom_num);

    return 0;
}

/*FUNCTIONS================================================================================================================*/
void read_input(){
    //user input for periodicity
    printf("periodicity along x: \n");
    scanf("%d", &ux);
    printf("periodicity along y: \n");
    scanf("%d", &uy);
    printf("periodicity along z: \n");
    scanf("%d", &uz);
    printf("input lattice constant\n");
    scanf("%lf", &a_constant);
}

void file_writing(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "LAB1_%s_neighb_list.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%.0lf \t %.0lf \t %lf \n",arr[i][0], arr[i][1] , arr[i][2]);
    }
    printf("neighbour list of %s has been printed!\n",lattice_structure);
    fclose(filepointer);
}

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

void fn_cross_product(int idx1 , int idx2 , int idx3){ //b_(idx1 + 1) = a_(idx2 + 1) x a_(idx3 + 1)
        b[idx1][0] = a[idx2][1]*a[idx3][2] - a[idx2][2]*a[idx3][1];
        b[idx1][1] = a[idx2][2]*a[idx3][0] - a[idx2][0]*a[idx3][2];
        b[idx1][2] = a[idx2][0]*a[idx3][1] - a[idx2][1]*a[idx3][0];
    return;
}

double fn_2pi_triple_product(){
    volume = (a[0][0] * b[0][0] + a[0][1] * b[0][1] + a[0][2] * b[0][2]);
    double answer = (2 * M_PI)/volume;
    return answer;
}

/*ATOM COORDS GENERATION==============================================================================================*/
void FN_atom_coords(){
    atom_num = ux * uy * uz;
    bcc_atom_num = atom_num * 2;
    fcc_atom_num = atom_num * 4;
    diamond_atom_num = atom_num * 8;

    //initialise 2D array for atom coords (initialise a atom_num x 3 array)
    sc_atom_coords = (double **)malloc(atom_num* sizeof(double *));
    for (int i = 0; i < atom_num; i++){
        sc_atom_coords[i] = (double *)malloc(3 * sizeof(double));
    }

    bcc_atom_coords = (double **)malloc(bcc_atom_num* sizeof(double *));
    for (int i = 0; i < bcc_atom_num; i++){
        bcc_atom_coords[i] = (double *)malloc(3 * sizeof(double));
    }

    fcc_atom_coords = (double **)malloc(fcc_atom_num* sizeof(double *));
    for (int i = 0; i < fcc_atom_num; i++){
        fcc_atom_coords[i] = (double *)malloc(3 * sizeof(double));
    }

    diamond_atom_coords = (double **)malloc(diamond_atom_num* sizeof(double *));
    for (int i = 0; i < diamond_atom_num; i++){
        diamond_atom_coords[i] = (double *)malloc(3 * sizeof(double));
    }

    fn_simplecubic(sc_atom_coords);
    fn_bcc(bcc_atom_coords);
    fn_fcc(fcc_atom_coords);
    fn_diamond_fcc(diamond_atom_coords);
    return;
}

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

void fn_bcc(double **atom_coords){
    fn_simplecubic(atom_coords); //construct a simple cubic lattice
    fn_shift_coords(atom_coords, atom_num, 0.5); //adding an atom to each basis atom at (a/2 , a/2 , a/2) relative to the basis atom
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

void fn_diamond_fcc(double **atom_coords){
    fn_fcc(atom_coords);
    //an atom placed at (a/4 , a/4 , a/4) relative to each fcc atom
    fn_shift_coords(atom_coords, atom_num * 4, 0.25);
    return;
}

/*NEIGHBOUR LIST GENERATION=======================================================================================*/
void FN_generate_neighbour_list(char *name, double** coords, int size){
    double row_num = size * (size - 1);
    double **neighbour_list = (double **)malloc(row_num * sizeof(double *));
    for (int i = 0; i < row_num; i++){
        neighbour_list[i] = (double *)malloc(3 * sizeof(double));
    }

    double distance_cutoff = fn_distance_cutoff(name);
    int neighbour_num = fn_neighbour_list(neighbour_list, coords, size,distance_cutoff);
    file_writing(name, neighbour_list, neighbour_num);

    free(neighbour_list);
    free(coords);
}

double fn_pbc(double* t){
    double n[3];
    double difference_plus = 0 , difference_minus = 0;
    for(int i = 0 ; i < 3 ; i++){
        n[i] = fmod((t[0] * b[i][0] + t[1] * b[i][1] + t[2] * b[i][2])/(2 * M_PI) , 1); //generate fractional numbers
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

double fn_distance_cutoff(char *name){ //calculate the distance cutoff for each type of structure
    double distance_cutoff = 0;
    if(name == SC_STR){
        distance_cutoff = (a_constant) * 1.001;
    }
    if(name == BCC_STR){
        distance_cutoff = (a_constant * sqrt(3) * 0.5) * 1.001;
    }
    if(name == FCC_STR){
        distance_cutoff = (a_constant * sqrt(2) * 0.5) * 1.001;
    }
    if(name == DIAMOND_FCC_STR){
        distance_cutoff = (a_constant * sqrt(3) * 0.25) * 1.001;
    }
    return distance_cutoff;
}