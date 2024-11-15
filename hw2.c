#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double s[4][4][4];

double FN_trace(double** arr){
    double sum = 0;
    for(int i = 0 ; i < 4;  i++){
        sum += arr[i][i];
    }
    return sum;
}

double** FN_product(double** arr_a ,  double** arr_b){
    double arr_product[4][4];

    for(int i = 0; i< 4 ; i++){
        for(int j = 0 ; j < 4 ; j++){
            for(int k = 0 ; k < 4 ;k++){
                arr_product[i][j] += arr_a[i][k] * arr_b[k][j];
            }
        }
    }

    return arr_product;
}

double** FN_inverse(double*** s_arr , double* arr_c){
    double d[4][4];
    for(int i = 0; i< 4 ; i++){
        for(int j = 0 ; j < 4 ; j++){
            d[i][j] = -s_arr[3][i][j] / arr_c[0];
        }
    }
    return d;
}

double*** FN_FL_algorithm(double** arr_a , double* arr_c){
    double s[4][4][4];

    for(int i = 0 ; i < 4 ; i++){s[0][i][i] = 1;}
    for(int k = 1; k < 4 ; k++){
        s[k] = FN_product(arr_a,s[k-1]);

    }

    return s;
}

int main(){

 
    double a[4][4] = { {1, 2, 3, 4 }, {5, 6, 7 , 8} , {9 , 10 , 11 , 12} , {13 , 14 , 15 , 16}};
    double c[4];

    FN_FL_algorithm(a,c);
    double** d = FN_inverse(s,c);

    for (int i = 0; i < 4 ; i++){
        printf("%lf \t %lf \t %lf \t %lf \n",d[i][0], d[i][1] , d[i][2], d[i][3], d[i][4]);
    }


    return 0;
}