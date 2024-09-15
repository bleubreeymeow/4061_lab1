#include <stdio.h>

double A[3];
double B[3];
double crossproduct[3];

double fn_dotproduct(int length){
    double dotproduct;
    for(int i = 0 ; i < length ; i++){
    dotproduct += A[i]*B[i];
    }
    return dotproduct;
}

void fn_crossproduct(){
    crossproduct[1] = A[2]*B[3] - A[3]*B[2];
    crossproduct[2] = A[3]*B[1] - A[1]*B[3];
    crossproduct[3] = A[1]*B[2] - A[2]*B[1];

}

int main(){

    printf("3 numbers for A:");
    for (int i = 0 ; i < 3 ; i++){
        scanf("%lf", &A[i]);
    }
    
       printf("3 numbers for B:");
    for (int i = 0 ; i < 3 ; i++){
        scanf("%lf", &B[i]);
    }
    
    int length = 3;

    double dotproduct = fn_dotproduct(length);
    fn_crossproduct();
    
    printf("dot product : %lf \n cross product : \n %lf \n %lf \n %lf \n ",dotproduct,crossproduct[1],crossproduct[2],crossproduct[3]);
    return 0;
}