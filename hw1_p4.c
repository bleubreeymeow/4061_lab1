#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double function(double x){
return x*x*x*x*x*x*x*x;
}

double fn_trapezoid(double h,double a,double b,int N){
    double sum_y;
    for ( int i = 1 ; i < N ; i++){
        double x = (a + h*i);
        sum_y += 2*function(x) ;
    }

    double y_0 = function(a);
    double y_N = function(b);

    double answer = 0.5*h*(y_0 + sum_y + y_N);

    return answer;
}

double fn_randomsample(double a , double b, int M){
    int i;
    double x , sum_y, answer,y,shift;

    sum_y = 0.;
    for (i = 0; i < M ; i++){
        x = a + (b - a) * ((double)rand() / RAND_MAX);
        sum_y += function(x);
    }
    
    answer = ((b-a)/M)*sum_y;
    return answer;
}

int main(){
    double a = -1.;
    double b =1.;
    int M = 100000;
    int N = 1000;
    double true_value  = 2./9;

    double h = (double)(b - a)/N ;

    double integrated_trapezoid = fn_trapezoid(h,a,b,N);
    double integrated_random = fn_randomsample(a,b,M);


    printf("%lf\n%lf",integrated_trapezoid, integrated_random);

    return 0;
}