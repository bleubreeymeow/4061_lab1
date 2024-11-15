#include <stdio.h>

int main() {
    int myInteger = 42;              // An integer value
    long double myLongDouble = 0.0;        // Declare a long double variable

    // Implicit conversion from int to long double
    myLongDouble = myInteger;

    // Print the values
    printf("Integer: %d\n", myInteger);
    printf("Long Double: %.20Lf\n", myLongDouble);

    return 0;
}