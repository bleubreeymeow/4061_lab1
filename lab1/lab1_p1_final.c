#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SC_STR "sc"
#define BCC_STR "bcc"
#define FCC_STR "fcc"
#define DIAMOND_FCC_STR "diamond"
#define ELEMENT_STR "Si"

int ux; //periodicity along x 
int uy; //periodicity along y
int uz; // periodicity along z
int atom_num; //atom number for sc lattice
double a; //lattice constant

void file_writing(char *lattice_structure, double **arr, int num);
void shift_coords(double **arr, int count, double factor);
void fn_simplecubic(double **arr);
void fn_bcc(double **arr);
void fn_fcc(double **arr);
void fn_diamond_fcc(double **arr);

int main(){
    //user input for periodicity
    printf("input ux: \n");
    scanf("%d", &ux);
    printf("input uy: \n");
    scanf("%d", &uy);
    printf("input uz: \n");
    scanf("%d", &uz);
    printf("input a: \n");
    scanf("%lf", &a);

    atom_num = ux * uy * uz;

    //initialise 2D array for lattice. rows = sc atom number*8 for maximum number of atoms (diamond structure)
    double **arr = (double **)malloc(atom_num * 8 * sizeof(double *));
    for (int i = 0; i < atom_num * 8; i++){
        arr[i] = (double *)malloc(3 * sizeof(double));
    }

    fn_simplecubic(arr);
    file_writing(SC_STR, arr, atom_num);

    fn_bcc(arr);
    file_writing(BCC_STR, arr, atom_num * 2);

    fn_fcc(arr);
    file_writing(FCC_STR, arr, atom_num * 4);

    fn_diamond_fcc(arr);
    file_writing(DIAMOND_FCC_STR, arr, atom_num * 8);
    printf("The structure of sc, bcc, fcc and diamond crystalline structures has been outputted!");
    free(arr);

    return 0;
}

void file_writing(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "lab1_%s.xyz", lattice_structure); //create filename
    FILE *filepointer = fopen(lattice_filename, "w");
    fprintf(filepointer,"%d \n \n",(num));
    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%s \t %lf \t %lf \t %lf \n", ELEMENT_STR, arr[i][0], arr[i][1] , arr[i][2]);
    }
    fclose(filepointer);
}

//shifting coords, for translation
void shift_coords(double **arr, int count, double factor){
    for (int i = 0; i < count; i++) {
        for (int j = 0; j < 3; j++) {
            arr[i + count][j] = arr[i][j] + a * factor;
        }
    }
}

void fn_simplecubic(double **arr) {
    for (int i = 0; i < atom_num ; i++) {
        arr[i][0] = (i % ux) * a;
        arr[i][1] = (i / ux) % uy * a;
        arr[i][2] = (i / ux / uy) % uz * a;
    }
    return;
}

void fn_bcc(double **arr){
    fn_simplecubic(arr); //construct a simple cubic lattice
    shift_coords(arr, atom_num, 0.5); //adding an atom to each basis atom at (a/2 , a/2 , a/2) relative to the basis atom
    return;
}

void fn_fcc(double **arr){
    fn_simplecubic(arr); //construct a simple cubic lattice
    for (int i = atom_num; i < atom_num * 4; i += 3) {
        int index = (int)((i - atom_num) / 3); //retrieve the atom index in the simple cubic lattice
        //2nd basis atom at (a/2 , a/2 , 0) relative to 1st basis atom
        arr[i][0] = arr[index][0] + a * 0.5;
        arr[i][1] = arr[index][1] + a * 0.5;
        arr[i][2] = arr[index][2];
        //3rd basis atom at (a/2 , 0 , a/2) relative to 1st basis atom       
        arr[i + 1][0] = arr[index][0] + a * 0.5;
        arr[i + 1][1] = arr[index][1];
        arr[i + 1][2] = arr[index][2] + a * 0.5;
        //4th basis atom at (0 , a/2 , a/2) relative to 1st basis atom
        arr[i + 2][0] = arr[index][0];
        arr[i + 2][1] = arr[index][1] + a * 0.5;
        arr[i + 2][2] = arr[index][2] + a * 0.5;
    }
    return;
}

void fn_diamond_fcc(double **arr){
    fn_fcc(arr);
    //an atom placed at (a/4 , a/4 , a/4) relative to each fcc atom
    shift_coords(arr, atom_num * 4, 0.25);
    return;
}