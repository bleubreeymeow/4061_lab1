#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void fn_cross_product(double **a, double **b, int idx1 , int idx2 , int idx3){ //b_(idx1 + 1) = a_(idx2 + 1) x a_(idx3 + 1)
        b[idx1][0] = a[idx2][1]*a[idx3][2] - a[idx2][2]*a[idx3][1];
        b[idx1][1] = a[idx2][2]*a[idx3][0] - a[idx2][0]*a[idx3][2];
        b[idx1][2] = a[idx2][0]*a[idx3][1] - a[idx2][1]*a[idx3][0];
    return;
}

double fn_2pi_triple_product(double **a , double **b){
    double answer = (2 * M_PI)/(a[0][0] * b[0][0] + a[0][1] * b[0][1] + a[0][2] * b[0][2]);
    return answer;
}

void fn_reciprocal(double **a, double **b){
    fn_cross_product(a , b ,0 , 1 , 2); // b1 = a2 x a3
    double scalars = fn_2pi_triple_product(a , b); //a1 dot (a2 x a3)
    fn_cross_product(a , b ,1 , 2 , 0); // b2 = a3 x a1
    fn_cross_product(a , b ,2 , 0 , 1); // b3 = a1 x a2

    for (int i = 0 ; i < 3 ; i++){
        for(int j = 0 ; j < 3 ; j++){
            b[i][j] = b[i][j] * scalars;
        }
    }
    return;
}

int main(){
    //initialise 2D array for a1 , a2 , a3 (initialise a 3x3 array)
    double **a = (double **)malloc(3 * sizeof(double *));
    for (int i = 0; i < 3 ; i++){
        a[i] = (double *)malloc(3 * sizeof(double));
    }
    //initialise 2D array for b1 , b2 , b3 (initialise a 3x3 array)
    double **b = (double **)malloc(3 * sizeof(double *));
    for (int j = 0; j < 3 ; j++){
        b[j] = (double *)malloc(3 * sizeof(double));
    }

    printf("enter a1 vector a1 xhat , a1 yhat , a1 zhat: \n");
    scanf("%lf %lf %lf", &a[0][0] , &a[0][1] , &a[0][2]);
    printf("enter a2 vector a2 xhat , a2 yhat , a2 zhat: \n");
    scanf("%lf %lf %lf", &a[1][0] , &a[1][1] , &a[1][2]);
    printf("enter a3 vector a3 xhat , a3 yhat , a3 zhat: \n");
    scanf("%lf %lf %lf", &a[2][0] , &a[2][1] , &a[2][2]);

    fn_reciprocal(a,b); //calculate the reciprocal vectors

    for( int k = 0 ; k < 3 ; k++){
        printf("b%d vector = %lf \t %lf \t %lf \n",k, b[k][0],  b[k][1],b[k][2]);
    }
    return 0;
}