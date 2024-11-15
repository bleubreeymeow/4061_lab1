#include <stdio.h>

int main() {
    int myInteger = 42;              // An integer value
    long double myLongDouble = 1.0;
    
    printf("%llf\n",myLongDouble);        // Declare a long double variable

    // Implicit conversion from int to long double
    myLongDouble = (long double) myInteger;

    // Print the values
    printf("Integer: %d\n", myInteger);
    printf("Long Double: %.20Lf\n", myLongDouble);

    return 0;
}