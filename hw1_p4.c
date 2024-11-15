#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double function(double x){
return x*x*x*x*x*x*x*x;
}

long double fn_trapezoid(long double h,double a,double b,long double N){
    long double sum_y;
    for ( int i = 1 ; i < N ; i++){
        double x = (a + h*i);
        sum_y += 2*function(x) ;
    }
    double y_0 = function(a);
    double y_N = function(b);

    return 0.5*h*(y_0 + sum_y + y_N);
}

long double fn_randomsample(double a , double b, int M){
    int i;
    long double x , sum_y;

    sum_y = 0.;
    for (i = 0; i < M ; i++){
        x = a + (b - a) * ((double)rand() / RAND_MAX);
        sum_y += function(x);
    }
    return ((b-a)/M)*sum_y;
}

int main(){
    double a = -1.;
    double b =1.;
    long double N = 100;
    long double true_value  = 2./9;

    long double h = (b - a)/N ;

    long double integrated_trapezoid = fn_trapezoid(h,a,b,N);
    long double integrated_random = fn_randomsample(a,b,N);


    long double sd_random = (integrated_random - true_value) * (integrated_random - true_value);
    long double sd_trapezoid  = (integrated_trapezoid - true_value) * (integrated_trapezoid - true_value);

    printf("%Le \t%.3Lf \t%.3Lf \t%.3Le \t%.3Le \n",N, integrated_trapezoid, integrated_random, sd_random, sd_trapezoid);
    //printf("%lf\n%lf\n%lf\n",integrated_trapezoid, integrated_random);

    return 0;
}