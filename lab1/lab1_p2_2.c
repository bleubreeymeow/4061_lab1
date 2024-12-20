#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double volume;
double b[3][3]; //reciprocal unit cell vectors
double a[3][3]; //unit cell vectors
double n[3];
double t[3];

void fn_cross_product(int idx1 , int idx2 , int idx3);
double fn_2pi_triple_product();
void FN_reciprocal();
void FN_pbc();
void read_input();

int main(){
    read_input();
    FN_reciprocal(); //calculate the reciprocal vectors and the volume
    FN_pbc();

    printf("fractional coordinates = %lf \t %lf \t %lf \n", n[0], n[1], n[2]);
    printf("x2 y2 z2 = %lf \t %lf \t %lf \n", t[0], t[1], t[2]);

    return 0;
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

void FN_pbc(){
    double difference_plus = 0 , difference_minus = 0;
    for(int i = 0 ; i < 3 ; i++){
        n[i] = fmod((t[0] * b[i][0] + t[1] * b[i][1] + t[2] * b[i][2])/(2 * M_PI) , 1); //generate fractional numbers

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
    return ;
}

void read_input(){
    //user input a1 a2 a3
    printf("enter a1 vector a1 xhat , a1 yhat , a1 zhat: \n");
    scanf("%lf %lf %lf", &a[0][0] , &a[0][1] , &a[0][2]);
    printf("enter a2 vector a2 xhat , a2 yhat , a2 zhat: \n");
    scanf("%lf %lf %lf", &a[1][0] , &a[1][1] , &a[1][2]);
    printf("enter a3 vector a3 xhat , a3 yhat , a3 zhat: \n");
    scanf("%lf %lf %lf", &a[2][0] , &a[2][1] , &a[2][2]);
    printf("enter x1, y1 , z1: \n");
    scanf("%lf %lf %lf", &t[0] , &t[1] , &t[2]);
    return;
}

