#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double trace(double (*a)[4]) {
    double sum = 0;
    for (int i = 0; i < 4; ++i) sum += a[i][i];
    return sum;
}

void mm(double (*c)[4], double (*a)[4], double (*b)[4]) {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                c[i][j] += a[i][k] * b[k][j];
    return;
}

void fl(double (*s)[4][4], double (*a)[4], double *c) {
    for (int i = 0; i < 4; ++i) s[0][i][i] = 1;
    for (int k = 1; k < 4; ++k) {
        mm(s[k], a, s[k - 1]);

        c[4 - k] = -trace(s[k]) / k;
        for (int i = 0; i < 4; ++i)
            s[k][i][i] += c[4 - k];
    }

    double temp[4][4] = { };
    mm(temp, a, s[4 - 1]);
    c[0] = -trace(temp) / 4;
    return;
}

int main() {

    printf("inverted matrix:\n");
    double a[4][4] = { {1, 2, 3, 4 }, {0, 1, 0 , 2} , {0 , 0 , 1 , 3} , {0 , 0 , 0 , 1}};


    double c[4] = { };
    double d[4][4] = { };
    double s[4][4][4] = { };

    fl(s, a, c);

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            d[i][j] = -s[4-1][i][j] / c[0];

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j)
            printf("%lf\t", d[i][j]);

        printf("\n");
    }
}