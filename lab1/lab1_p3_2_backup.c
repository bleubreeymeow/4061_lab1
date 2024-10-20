#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define FCC_STR "fcc2_MgO"
#define COULOMB 1./(4.*M_PI*0.005526349406) //8.9875
#define ELEMENTARY_CHARGE 1//1.602E-19
#define A_Mg_Mg 423.989 //0
#define B_Mg_Mg 3.15
#define RHO_Mg_Mg 1./3.15 //0
#define C_Mg_Mg 1.049 //0

#define A_Mg_O  1216.35 //821.6
#define B_Mg_O 3.15
#define RHO_Mg_O 1./3.15 //0.3242
#define C_Mg_O 8.71 //0

#define A_O_O 3489.506//22764.0
#define B_O_O 3.15
#define RHO_O_O 1./3.15// 0.149
#define C_O_O  72.3708//27.88

int ux = 5, uy = 5 , uz = 5;
double a_constant = 5.578; //5.61;
double distance_cutoff = 5.578 * 5; //19.9; //maybe future caluclations for different cutoff distance from lattice constant
int atom_num;
double volume;
double b[3][3]; //reciprocal unit cell vectors
double a[3][3]; //unit cell vectors
double **sc_atom_coords = NULL;
double **sc_neighbour_list = NULL;

/*FUNCTIONS================================================================================================================*/

void file_writing(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "ACOULOMB_%s_neighb_list.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
       fprintf(filepointer,"%.0lf \t %.0lf \t %lf \t %lf \t %lf\n",arr[i][0], arr[i][1] , arr[i][2], arr[i][3], arr[i][4]);
    }
    printf("neighbour list of %s has been printed!\n",lattice_structure);
    fclose(filepointer);
}

void file_writing_coords(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "ACOULOMBLAB1_%s_atom_coords.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%lf \t %lf \t %lf \t %.0lf\n",arr[i][0], arr[i][1] , arr[i][2], arr[i][3]);
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
    int atom_type = 1;
    for (int i = 0; i < atom_num ; i++) {
        atom_coords[i][0] = (i % ux) * a_constant;
        atom_coords[i][1] = (i / ux) % uy * a_constant;
        atom_coords[i][2] = (i / ux / uy) % uz * a_constant;

        atom_type *= -1;
        if(ux % 2 == 0 && i % (ux * uy) == 0){
            atom_type *= -1;
        }
        if(ux % 2 == 0 && i % ux == 0){ //if perodicity of y is an even number and the loop has reached to the end, then the atom in the next y 
            atom_type *= -1;
        }
        atom_coords[i][3] = atom_type; //all the 'simple cubic' lattice points are Mg
  
    }
    return;
}

void FN_atom_coords(){
    atom_num = ux * uy * uz;
    sc_atom_coords = (double **)malloc(atom_num* sizeof(double *));
    for (int i = 0; i < atom_num; i++){
        sc_atom_coords[i] = (double *)malloc(4 * sizeof(double));
    }
    fn_simplecubic(sc_atom_coords);
    return;
}

/*NEIGHBOUR LIST GENERATION=======================================================================================*/
double fn_pbc(double* t){
    double n[3];
    double difference_plus = 0 , difference_minus = 0;
    for(int i = 0 ; i < 3 ; i++){
        n[i] = fmod((t[0] * b[i][0] + t[1] * b[i][1] + t[2] * b[i][2]) , 1); //generate fractional numbers
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
    int interaction_type = 0; //if both are Mg atoms, interaction = -2, if one is Mg , one is O, interaction = 0 , if both oxygen, interaction = 2

    for(int i = 0 ; i < size ; i++){
        for(int j = i + 1 ; j < size ; j++){ //j = i + 1 to avoid double counting (1 2 , 2 1), i + 1 for skipping (0 0) , (1 1)
            t[0] = source[j][0] - source[i][0]; //dx
            t[1] = source[j][1] - source[i][1]; //dy
            t[2] = source[j][2] - source[i][2]; //dz
            interaction_type = source[j][3] + source[i][3];

            distance = fn_pbc(t);

            if(distance <= distance_cutoff){ //evaulate whether atom i and atom j are nearest neighbours
                target[nearest_neighbour_num][0] = i;
                target[nearest_neighbour_num][1] = j;
                target[nearest_neighbour_num][2] = distance;
                target[nearest_neighbour_num][3] = interaction_type;
                nearest_neighbour_num = nearest_neighbour_num + 1;
            }
        }
    }
    return nearest_neighbour_num;
}


double fn_coulomb_buck_potential(double** neighbour_list, int size){
    double potential = 0;
    double temp_potential = 0;
    double r_6 = 0;
    double r = 0; 

    for(int i = 0 ; i < size ; i++){
        //r_6 = neighbour_list[i][2]* neighbour_list[i][2] *  neighbour_list[i][2]* neighbour_list[i][2] *  neighbour_list[i][2]* neighbour_list[i][2];
        r = neighbour_list[i][2];
        r_6 = r * r * r * r * r * r;
        if(neighbour_list[i][3] == 2){ //Mg-Mg interaction
            temp_potential = COULOMB / r + A_Mg_Mg * exp(-r * B_Mg_Mg) - C_Mg_Mg / r_6 ;
            neighbour_list[i][4] = temp_potential;
        }
        if(neighbour_list[i][3] == 0){ //Mg-O interaction
            temp_potential = (-COULOMB  / (r))+ (A_Mg_O * exp(-r * B_Mg_O)) - (C_Mg_O / r_6);
            neighbour_list[i][4] = temp_potential;
        }
        if(neighbour_list[i][3] == -2){ //O-O interaction
            temp_potential = ( COULOMB / r) + A_O_O * exp(-r * B_O_O) - (C_O_O / r_6);
            neighbour_list[i][4] = temp_potential;
        }
        printf("%lf\n",temp_potential);
        potential += temp_potential;

    }
    return potential;
}

void FN_generate_neighbour_list(char *name, double** coords, int size){
    double row_num = size * (size - 1);
    double **neighbour_list = (double **)malloc(row_num * sizeof(double *));
    for (int i = 0; i < row_num; i++){
        neighbour_list[i] = (double *)malloc(5 * sizeof(double));
    }

    int neighbour_num = fn_neighbour_list(neighbour_list, coords, size);
    double coulomb_buck_potential = fn_coulomb_buck_potential(neighbour_list, neighbour_num);
    printf("\n %s coulomb buckingham potential(eV): %e\n", name, coulomb_buck_potential/(atom_num));
    file_writing(name, neighbour_list, neighbour_num);

    free(neighbour_list);
    free(coords);
}


/*MAIN==============================================================================================================*/
int main(){
    FN_a_vect(); //initialise a1 a2 a3 cell vectors
    FN_reciprocal(); //calculate the reciprocal vectors and the volume
    FN_atom_coords(); //generate fcc coordinates
    file_writing_coords(FCC_STR,sc_atom_coords,atom_num);
    //generate fcc neighbour list
    FN_generate_neighbour_list(FCC_STR, sc_atom_coords, atom_num);

    return 0;
}