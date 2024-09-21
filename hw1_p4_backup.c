#include<stdio.h>
#include<stdlib.h>
#include<math.h>

long double function(long double x);
long double fn_trapezoid(long double h, long double a, long double b, long double N);
long double fn_randomsample(double a ,double b, long double M);

int main(){
    double a = -1.;
    double b = 1;
    long double N = (long double)1;
    long double true_value  = (long double)2/9;
    int iteration_num = 10;
    long double h = (long double)0;

    for (int i = 0; i < iteration_num; i++){
        N *= 10;
        h = (long double)(b - a)/N ;

        long double integrated_trapezoid = fn_trapezoid(h,a,b,N);
        long double integrated_random = fn_randomsample(a,b,N);

        long double sd_random = (integrated_random - true_value) * (integrated_random - true_value);
        long double sd_trapezoid  = (integrated_trapezoid - true_value) * (integrated_trapezoid - true_value);

        printf("%.2Le \t%.5Lf \t%.5Lf \t%.3Le \t%.3Le \n",N, integrated_trapezoid, integrated_random, sd_trapezoid, sd_random);
    }
    return 0;
}

long double function(long double x){
return x*x*x*x*x*x*x*x;
}

long double fn_trapezoid(long double h, long double a, long double b, long double N){
    long double sum_y = (long double)0;

    for ( int i = 1 ; i < N ; i++){
        long double x = (a + h*(long double)i);
        sum_y += 2*function(x) ;
    }

    long double y_0 = function(a);
    long double y_N = function(b);

    return 0.5*h*(y_0 + sum_y + y_N);
}

long double fn_randomsample(double a , double b, long double M){
    int i;
    long double x = (long double)0;
    long double sum_y = (long double)0;

    for (i = 0; i < M ; i++){
        x = a + (b - a) * ((long double)rand() / RAND_MAX);
        sum_y += function(x);
    }
    return ((b-a)/M)*sum_y;
}