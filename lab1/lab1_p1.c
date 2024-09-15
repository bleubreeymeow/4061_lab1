#include<stdio.h>
#include <stdlib.h>
#include <string.h>

#define SC_STR "sc"
#define BCC_STR "bcc"
#define FCC_STR "fcc"
#define DIAMOND_FCC_STR "diamond"
#define ELEMENT_STR "Si"

int ux; //periodicity
int uy;
int uz;
int atom_num;
double a;

void file_writing(char *lattice_structure, double **arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "lab1_%s.xyz", lattice_structure);

    FILE *filepointer = fopen(lattice_filename, "w");
    
    fprintf(filepointer,"%d \n \n",(num));

    for (int i = 0; i < num ; i++){
        fprintf(filepointer,"%s \t %lf \t %lf \t %lf \n", ELEMENT_STR, arr[i][0], arr[i][1] , arr[i][2]);
        //MATLAB fprintf(filepointer,"%lf \t %lf \t %lf \n",arr[i][0], arr[i][1] , arr[i][2]);
    }

    fclose(filepointer);
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
    for (int i = 0; i < atom_num * 2 ; i = i + 2) {
        //1st basis atom
        arr[i][0] = (i % ux) * a;
        arr[i][1] = (i / ux) % uy * a;
        arr[i][2] = (i / ux / uy) % uz * a;
        //2nd basis atom at (a/2, a/2, a/2) relative to 1st basis atom
        arr[i + 1][0] = arr[i][0] + a * 0.5;
        arr[i + 1][1] = arr[i][1] + a * 0.5;
        arr[i + 1][2] = arr[i][2] + a * 0.5;
    }
    return;
}

void fn_fcc(double **arr){
    for (int i = 0; i < atom_num * 4 ; i = i + 4) {
        //1st basis atom
        arr[i][0] = (i % ux) * a;
        arr[i][1] = (i / ux) % uy * a;
        arr[i][2] = (i / ux / uy) % uz * a;
        //2nd basis atom at (a/2 , a/2 , 0) relative to 1st basis atom
        arr[i + 1][0] = arr[i][0] + a * 0.5;
        arr[i + 1][1] = arr[i][1] + a * 0.5;
        arr[i + 1][2] = arr[i][2];
        //3rd basis atom at (a/2 , 0 , a/2) relative to 1st basis atom       
        arr[i + 2][0] = arr[i][0] + a * 0.5;
        arr[i + 2][1] = arr[i][1];
        arr[i + 2][2] = arr[i][2] + a * 0.5;
        //4th basis atom at (0 , a/2 , a/2) relative to 1st basis atom
        arr[i + 3][0] = arr[i][0];
        arr[i + 3][1] = arr[i][1] + a * 0.5;
        arr[i + 3][2] = arr[i][2] + a * 0.5;
    }
    return;
}

void fn_diamond_fcc(double **arr){
    for (int i = 0; i < atom_num * 8 ; i = i + 8) {
        //1st basis atom
        arr[i][0] = (i % ux) * a;
        arr[i][1] = (i / ux) % uy * a;
        arr[i][2] = (i / ux / uy) % uz * a;
        //(a/4 , a/4 , a/4) relative to 1st fcc basis atom
        arr[i + 1][0] = arr[i][0] + a * 0.25;
        arr[i + 1][1] = arr[i][1] + a * 0.25;
        arr[i + 1][2] = arr[i][2] + a * 0.25;

        //================
        //2nd basis atom at (a/2, a/2, 0) from 1st basis atom
        arr[i + 2][0] = arr[i][0] + a * 0.5;
        arr[i + 2][1] = arr[i][1] + a * 0.5;
        arr[i + 2][2] = arr[i][2];
        //(a/4 , a/4 , a/4) relative to 2nd fcc basis atom
        arr[i + 3][0] = arr[i + 2][0] + a * 0.25;
        arr[i + 3][1] = arr[i + 2][1] + a * 0.25;
        arr[i + 3][2] = arr[i + 2][2] + a * 0.25;

        //===================
        //3rd basis atom at (a/2, 0, a/2) from 1st basis atom
        arr[i + 4][0] = arr[i][0] + a * 0.5;
        arr[i + 4][1] = arr[i][1];
        arr[i + 4][2] = arr[i][2] + a * 0.5;
        //(a/4 , a/4 , a/4) relative to 3rd fcc basis atom
        arr[i + 5][0] = arr[i + 4][0] + a * 0.25;
        arr[i + 5][1] = arr[i + 4][1] + a * 0.25;
        arr[i + 5][2] = arr[i + 4][2] + a * 0.25;

        //====================
        //4th basis atom at (0, a/2, a/2) from 1st basis atom
        arr[i + 6][0] = arr[i][0];
        arr[i + 6][1] = arr[i][1] + a * 0.5;
        arr[i + 6][2] = arr[i][2] + a * 0.5;
        //(a/4 , a/4 , a/4) relative to 4th fcc basis atom
        arr[i + 7][0] = arr[i + 6][0] + a * 0.25;
        arr[i + 7][1] = arr[i + 6][1] + a * 0.25;
        arr[i + 7][2] = arr[i + 6][2] + a * 0.25;
    }
    return;
}

int main(){

    printf("input ux: \n");
    scanf("%d", &ux);

    printf("input uy: \n");
    scanf("%d", &uy);

    printf("input uz: \n");
    scanf("%d", &uz);

    printf("input a: \n");
    scanf("%lf", &a);

    atom_num = ux * uy * uz;

    //initialise 2D array for sc lattice
    double **sc_arr = (double **)malloc(atom_num * sizeof(double *));
    for (int i = 0; i < atom_num ; i++){
        sc_arr[i] = (double *)malloc(3 * sizeof(double));
    }

    //initialise 2D array for bcc lattice
    double **bcc_arr = (double **)malloc(atom_num * 2 * sizeof(double *));
    for (int j = 0; j < atom_num * 2 ; j++){
        bcc_arr[j] = (double *)malloc(3 * sizeof(double));
    }

    //initialise 2D array for fcc lattice
    double **fcc_arr = (double **)malloc(atom_num * 4 * sizeof(double *));
    for (int k = 0; k < atom_num * 4 ; k++){
        fcc_arr[k] = (double *)malloc(3 * sizeof(double));
    }

    //initialise 2D array for diamond lattice
    double **diamond_fcc_arr = (double **)malloc(atom_num * 8 * sizeof(double *));
    for (int m = 0; m < atom_num * 8 ; m++){
        diamond_fcc_arr[m] = (double *)malloc(3 * sizeof(double));
    }

    fn_simplecubic(sc_arr);
    file_writing(SC_STR, sc_arr,atom_num);
    fn_bcc(bcc_arr);
    file_writing(BCC_STR , bcc_arr,atom_num * 2);
    fn_fcc(fcc_arr);
    file_writing(FCC_STR , fcc_arr, atom_num * 4);
    fn_diamond_fcc(diamond_fcc_arr);
    file_writing(DIAMOND_FCC_STR ,diamond_fcc_arr, atom_num * 8);

    free(sc_arr);
    free(bcc_arr);
    free(fcc_arr);
    free(diamond_fcc_arr);

    return 0;
}