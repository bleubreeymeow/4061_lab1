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
        printf("n = %lf\n",n[i]);

    }


    for(int i = 0 ; i < 3 ; i++){
        t[i] = n[0] * a[0][i] + n[1] * a[1][i] + n[2] * a[2][i];
        printf("t = %lf\n",t[i]);
    }

    double distance = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
    return distance;

}



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

   distance_cutoff = 1.26;

 
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


    double t[3] = {3.75, 0, 0};

    FN_reciprocal(); //calculate the reciprocal vectors and the volume
    double distance = fn_pbc(t);

    printf("%lf\n",distance);


    for( int k = 0 ; k < 3 ; k++){
        printf("b%d vector = %lf \t %lf \t %lf \n",k, b[k][0],  b[k][1],b[k][2]);
    }
    printf("volume = %lf \n" , volume);


    
    return 0;
}
