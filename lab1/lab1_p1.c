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

void file_writing(char *lattice_structure, double *arr, int num){
    char lattice_filename[50];
    snprintf(lattice_filename, sizeof(lattice_filename), "lab1_%s.txt", lattice_structure);

    FILE *filepointer = fopen(lattice_filename, "w");
    
    //fprintf(filepointer,"%d \n \n",(num));

    for (int i = 0; i < num ; i++){
        //fprintf(filepointer,"%s \t %lf \t %lf \t %lf \n", ELEMENT_STR, arr[3 * i], arr[3 * i + 1] , arr[3 * i + 2]);
        fprintf(filepointer,"%lf \t %lf \t %lf \n", arr[3 * i], arr[3 * i + 1] , arr[3 * i + 2]);
    }

    fclose(filepointer);
}

void fn_simplecubic(double *arr) {
    for (int i = 0; i < atom_num ; i++) {
        arr[3 * i] = (i % ux) * a;
        arr[3 * i + 1] = (i / ux) % uy * a;
        arr[3 * i + 2] = (i / ux / uy) % uz * a;

    }
    return;
}

void fn_bcc(double *arr){
    for (int i = 0; i < atom_num ; i++) {
        arr[6 * i] = (i % (ux)) * (a) ;
        arr[6 * i + 1] = (i / ux) % uy * a;
        arr[6 * i + 2] = (i / ux / uy) % uz * a;

        arr[6 * i + 3] = arr[6 * i] + a * 0.5;
        arr[6 * i + 4] = arr[6 * i + 1] + a * 0.5;
        arr[6 * i + 5] = arr[6 * i + 2] + a * 0.5;
    }
    return;
}

void fn_fcc(double *arr){
    for (int i = 0; i < atom_num ; i++) {
        arr[12 * i] = (i % (ux)) * (a) ;
        arr[12 * i + 1] = (i / ux) % uy * a;
        arr[12 * i + 2] = (i / ux / uy) % uz * a;

        arr[12 * i + 3] = arr[12 * i] + a * 0.5;
        arr[12 * i + 4] = arr[12 * i + 1] + a * 0.5;
        arr[12 * i + 5] = arr[12 * i + 2];

        arr[12 * i + 6] = arr[12 * i] + a * 0.5;
        arr[12 * i + 7] = arr[12 * i + 1];
        arr[12 * i + 8] = arr[12 * i + 2] + a * 0.5;

        arr[12 * i + 9] = arr[12 * i];
        arr[12 * i + 10] = arr[12 * i + 1] + a * 0.5;
        arr[12 * i + 11] = arr[12 * i + 2] + a * 0.5;
    }
    return;
}

void fn_diamond_fcc(double *arr){
    for (int i = 0; i < atom_num ; i++) {
        //1st basis atom
        arr[24 * i] = (i % (ux)) * (a) ;
        arr[24 * i + 1] = (i / ux) % uy * a;
        arr[24 * i + 2] = (i / ux / uy) % uz * a;
        //(1/4 , 1/4 , 1/4) relative to 1st fcc basis atom
        arr[24 * i + 12] = arr[12 * i] + a * 0.25;
        arr[24 * i + 13] = arr[12 * i + 1] + a * 0.25;
        arr[24 * i + 14] = arr[12 * i + 2] + a * 0.25;

        //================
        //2nd basis atom at (1/2,1/2,0) from 1st basis atom
        arr[24 * i + 3] = arr[12 * i] + a * 0.5;
        arr[24 * i + 4] = arr[12 * i + 1] + a * 0.5;
        arr[24 * i + 5] = arr[12 * i + 2];
        //(1/4 , 1/4 , 1/4) relative to 2nd fcc basis atom
        arr[24 * i + 15] = arr[12 * i + 3] + a * 0.25;
        arr[24 * i + 16] = arr[12 * i + 4] + a * 0.25;
        arr[24 * i + 17] = arr[12 * i + 5] + a * 0.25;

        //===================
        //3rd basis atom at (1/2,0,1/2) from 1st basis atom
        arr[24 * i + 6] = arr[12 * i] + a * 0.5;
        arr[24 * i + 7] = arr[12 * i + 1];
        arr[24 * i + 8] = arr[12 * i + 2] + a * 0.5;
        //(1/4 , 1/4 , 1/4) relative to 3rd fcc basis atom
        arr[24 * i + 18] = arr[12 * i + 6] + a * 0.25;
        arr[24 * i + 19] = arr[12 * i + 7] + a * 0.25;
        arr[24 * i + 20] = arr[12 * i + 8] + a * 0.25;

        //====================
        //4th basis atom at (0,1/2,1/2) from 1st basis atom
        arr[24 * i + 9] = arr[12 * i];
        arr[24 * i + 10] = arr[12 * i + 1] + a * 0.5;
        arr[24 * i + 11] = arr[12 * i + 2] + a * 0.5;
        //(1/4 , 1/4 , 1/4) relative to 4th fcc basis atom
        arr[24 * i + 21] = arr[12 * i + 9] + a * 0.25;
        arr[24 * i + 22] = arr[12 * i + 10] + a * 0.25;
        arr[24 * i + 23] = arr[12 * i + 11] + a * 0.25;


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

    double *sc_arr = (double *)malloc(atom_num * 3 * sizeof(double));
    double *bcc_arr = (double *)malloc(atom_num * 6 * sizeof(double));
    double *fcc_arr = (double *)malloc(atom_num * 12 * sizeof(double));
    double *diamond_fcc_arr = (double *)malloc(atom_num * 24 * sizeof(double));

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