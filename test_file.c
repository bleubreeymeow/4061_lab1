#include<stdio.h>
int main(){
    int number1, number2, sum;

    printf("print two integers:\n");
    scanf("%d %d", &number1 ,&number2);

    sum = number1 + number2;

    printf("sum is %d \n", sum);

    return 0 ;
}