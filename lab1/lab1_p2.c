#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ux , uy , uz;
double AX , AY , AZ;
double a_constant;
double distance_cutoff;
double structure;
int atom_num;
int bcc_atom_num;
int fcc_atom_num;
int diamond_atom_num;
double volume;
double b[3][3];
double a[3][3] = { {5 , 0 , 0}, {0 , 5 , 0} , {0 , 0, 5}};
double **sc_atom_coords = NULL;
double **bcc_atom_coords = NULL;
double **fcc_atom_coords = NULL;
double **diamond_atom_coords = NULL;

double **sc_neighbour_list = NULL;
double **bcc_neighbour_list = NULL;
double **fcc_neighbour_list = NULL;
double **diamond_neighbour_list = NULL;

#define SC_STR "sc2"
#define BCC_STR "bcc2"
#define FCC_STR "fcc2"
#define DIAMOND_FCC_STR "diamond2"
#define ELEMENT_STR "Si2"


void file_writing_atom_coords(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "lab1_%s.txt", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    //fprintf(filepointer,"%d \n \n",(num));
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%lf \t %lf \t %lf \n",arr[i][0], arr[i][1] , arr[i][2]);
    }
    fclose(filepointer);
}

/*b1 b2 b3 VECTOR GENERATION=============================================================================================*/

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

void fn_simplecubic(double **atom_coords) {
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

void FN_atom_coords(){
    atom_num = ux * uy * uz;
    bcc_atom_num = atom_num * 2;
    fcc_atom_num = atom_num * 4;
    diamond_atom_num = atom_num * 8;

    //initialise 2D array for atom coords (initialise a (atom_num * structure_factor) x 3 array)
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


/*NEIGHBOUR LIST GENERATION=======================================================================================*/

double fn_pbc(double* t){
    double n[3];
    double difference_plus = 0 , difference_minus = 0;
    for(int i = 0 ; i < 3 ; i++){
        
        n[i] = fmod((t[0] * b[i][0] + t[1] * b[i][1] + t[2] * b[i][2])/(2 * M_PI) , 1);

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
        t[i] = n[0] * a[0][i] + n[1] * a[1][i] + n[2] * a[2][i];
    }
    double distance = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
    return distance;
}

void FN_neighbour_list(){
    double distance = 0;
    int k = 0;
    double t[3];
    double row_num = (atom_num * atom_num) - atom_num;
    sc_neighbour_list = (double **)malloc(row_num * sizeof(double *));
    for (int i = 0; i < row_num; i++){
        sc_neighbour_list[i] = (double *)malloc(3 * sizeof(double));
    }

    for(int i = 0 ; i < atom_num ; i++){
        for(int j = 0 ; j < atom_num ; j++){
            if(i == j){continue;}
                t[0] = sc_atom_coords[j][0] - sc_atom_coords[i][0]; //dx
                t[1] = sc_atom_coords[j][1] - sc_atom_coords[i][1]; //dy
                t[2] = sc_atom_coords[j][2] - sc_atom_coords[i][2]; //dz

                distance = fn_pbc(t);
                
               if(distance < distance_cutoff){
                    
                    sc_neighbour_list[k][0] = i;
                    sc_neighbour_list[k][1] = j;
                    sc_neighbour_list[k][2] = distance;

                    k = k + 1;
              
                  }
        }

    }
   
    return;
}

/*MAIN==============================================================================================================*/

int main(){

    /*printf("enter a1 vector a1 xhat , a1 yhat , a1 zhat: \n");
    scanf("%lf %lf %lf", &a[0][0] , &a[0][1] , &a[0][2]);
    printf("enter a2 vector a2 xhat , a2 yhat , a2 zhat: \n");
    scanf("%lf %lf %lf", &a[1][0] , &a[1][1] , &a[1][2]);
    printf("enter a3 vector a3 xhat , a3 yhat , a3 zhat: \n");
    scanf("%lf %lf %lf", &a[2][0] , &a[2][1] , &a[2][2]);

    */

   ux = 4;
   uy = 4;
   uz = 4;
   a_constant = 1.25;
   distance_cutoff = 1.30;
 
    //user input for periodicity
    /*printf("input number of atoms along x in the unit cell: \n");
    scanf("%d", &ux);
    printf("input number of atoms along y in the unit cell: \n");
    scanf("%d", &uy);
    printf("input number of atoms along z in the unit cell: \n");
    scanf("%d", &uz);
    printf("input distance between basis atoms: \n");
    scanf("%lf", &a_constant);
    printf("input unit cell length along x: \n");
    scanf("%lf", &AX);
    printf("input unit cell length along y: \n");
    scanf("%lf", &AY);
    printf("input unit cell length along z: \n");
    scanf("%lf", &AZ);
    printf("input distance cutoff for neighbour list:\n");
    scanf("%lf", &distance_cutoff);
*/

    FN_reciprocal(); //calculate the reciprocal vectors and the volume
    FN_atom_coords(); //generate an array with atom coordinates

    file_writing_atom_coords(SC_STR, sc_atom_coords, atom_num);

    FN_neighbour_list();

    printf("neighbor sc \n");
    for( int ii = 0 ; ii < atom_num*6; ii++){
        printf("%lf \t %lf \t %lf \n", sc_neighbour_list[ii][0],  sc_neighbour_list[ii][1],sc_neighbour_list[ii][2]);
    }

    for( int k = 0 ; k < 3 ; k++){
        printf("b%d vector = %lf \t %lf \t %lf \n",k, b[k][0],  b[k][1],b[k][2]);
    }
    printf("volume = %lf \n" , volume);

    free(sc_atom_coords);
    free(bcc_atom_coords);
    free(fcc_atom_coords);
    free(diamond_atom_coords);
    
    return 0;
}
