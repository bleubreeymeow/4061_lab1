#include<stdio.h>
#include <stdlib.h>
#include <string.h>

int *atom_num;
double *a;

double *xarray_sc;
double *yarray_sc;
double *zarray_sc;

double *xarray_bcc;
double *yarray_bcc;
double *zarray_bcc;

double *xarray_fcc;
double *yarray_fcc;
double *zarray_fcc;

#define element = "Si"

#define SC = "sc"
#define BCC = "bcc"
#define FCC = "fcc"

void file_writing(char *lattice_structure, double xarray[] , double yarray[] , double zarray[]){
    char lattice_filename[25];
    snprintf(lattice_filename, sizeof(lattice_filename), "lab1_%s.xyz", lattice_structure);

    FILE *filepointer = fopen(lattice_filename, "w");
    
    
    fprintf(filepointer,"%d \n \n",(*atom_num));
    //fprintf(filepointer,"\n");
    for (int i = 0; i < *atom_num ; i++){
        fprintf(filepointer,"%s \t %lf \t %lf \t %lf \n", element, xarray[i], yarray[i] , zarray[i]);
    }
    fclose(filepointer);
}

void fn_simplecubic(){

    //allocate memory for arrays
    xarray_sc = (double *)malloc(*atom_num * sizeof(double));
    yarray_sc = (double *)malloc(*atom_num * sizeof(double));
    zarray_sc = (double *)malloc(*atom_num * sizeof(double));

    for (int i = 0; i < (*atom_num) ; i++) {
        xarray_sc[i] = (i%(ux))*(*a);
        yarray_sc[i] = (i/(ux))%(uy)*(*a);
        zarray_sc[i] = (i/(ux)/(uy))%(uz)*(*a);
        printf("%lf \t %lf \t %lf \n",xarray_sc[i],yarray_sc[i],zarray_sc[i]);
    }
    
    return;
}

/*void fn_bcc(){
    int bcc_atom_num = *atom_num*2;
    //allocate memory for arrays
    xarray_bcc = (double *)malloc(bcc_atom_num * sizeof(double));
    yarray_bcc = (double *)malloc(bcc_atom_num * sizeof(double));
    zarray_bcc = (double *)malloc(bcc_atom_num * sizeof(double));


    for (int i = 0; i < bcc_atom_num ; i++) {
        for(int j = 0; j < bcc_atom_num ; j++){
            for(int k = 0; k < bcc_atom_num; k++){
                xarray_bcc[i] = ((i-j+k)%(*ux))*(*a)*0.5;
                yarray_bcc[i] = ((i+j-k)/(*ux))%(*uy)*(*a)*0.5;
                zarray_bcc[i] = ((-i+j+k)/(*ux)/(*uy))%(*uz)*(*a)*0.5;
            }
        } 
    }
    
    return;
}
*/


int main(){
    int ux, uy, uz;
    double a_main;
    printf("input ux: \n");
    scanf("%d", &ux);
    printf("input uy: \n");
    scanf("%d", &uy);
    printf("input uz: \n");
    scanf("%d", &uz);
    printf("input a: \n");
    scanf("%lf", &a_main);

    a = &a_main;

    int mainatom_num = ux*uy*uz;
    atom_num = &mainatom_num;

    fn_simplecubic();
    file_writing(SC, xarray_sc,yarray_sc , zarray_sc);

    //fn_bcc();

    //file_writing(bcc,xarray_bcc,yarray_bcc , zarray_bcc);



    //free(xarray_bcc);
    //free(yarray_bcc);
    //free(zarray_bcc);


    free(xarray_sc);
    free(yarray_sc);
    free(zarray_sc);
    return 0;
}