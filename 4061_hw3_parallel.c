#include <math.h>
#include <stdio.h>

#define N 100
#define J 2

#define G 9.80
#define ANGLE 42.5 * M_PI / 180
#define SPEED 67.0
#define MASS 250.0
#define AREA 0.93
#define DENSITY 1.2
#define K AREA * DENSITY / (2 * MASS)
#define DT 2 * SPEED * sin(ANGLE) / (G * N)
#define D DT * DT / 2

int main(){
    double x[N + 1] = { 0 };
    double y[N + 1] = { 0 };
    double vx[N + 1] = { 0 };
    double vy[N + 1] = { 0 };
    double ax[N + 1] = { 0 };
    double ay[N + 1] = { 0 };

    vx[0] = SPEED * cos(ANGLE);
    vy[0] = SPEED * sin(ANGLE);

    double v = sqrt(vx[0] * vx[0] + vy[0] * vy[0]);

    ax[0] = -K * pow(v,3./4) * vx[0];
    ay[0] = -G - K * pow(v,3./4) * vy[0];

    double p = vx[0] * ax[0] + vy[0] * ay[0];

    x[1] = x[0] + DT * vx[0] + D * ax[0];
    y[1] = y[0] + DT * vy[0] + D * ay[0];

    vx[1] = vx[0] + DT * ax[0] - D * K * (pow(v,3./4) * ax[0] + p * (3./4) * vx[0] / pow(v,5./4));
    vy[1] = vy[0] + DT * ay[0] - D * K * (pow(v,3./4) * ay[0] + p * (3./4) * vy[0] / pow(v,5./4));

    v = sqrt(vx[1] * vx[1] + vy[1] * vy[1]);

    ax[1] = -K * pow(v,3./4) * vx[1];
    ay[1] = -G - K * pow(v,3./4) * vy[1];

    //Calculate other position and velocity recursively

    double d2 = 2 * DT;
    double d3 = DT / 3;

    for (int i = 0; i < N - 1; ++i) {
        //Predict the next position and velocity

        x[i + 2] = x[i] + d2 * vx[i + 1];
        y[i + 2] = y[i] + d2 * vy[i + 1];

        vx[i + 2] = vx[i] + d2 * ax[i + 1];
        vy[i + 2] = vy[i] + d2 * ay[i + 1];

        v = sqrt(vx[i + 2] * vx[i + 2] + vy[i + 2] * vy[i + 2]);
        
        ax[i + 2] = -K * pow(v,3./4) * vx[i + 2];
        ay[i + 2] = -G - K * pow(v,3./4) * vy[i + 2];

        //Correct the new position and velocity

        x[i + 2] = x[i] + d3 * (vx[i + 2] + 4 * vx[i + 1] + vx[i]);
        y[i + 2] = y[i] + d3 * (vy[i + 2] + 4 * vy[i + 1] + vy[i]);

        vx[i + 2] = vx[i] + d3 * (ax[i + 2] + 4 * ax[i + 1] + ax[i]);
        vy[i + 2] = vy[i] + d3 * (ay[i + 2] + 4 * ay[i + 1] + ay[i]);
    }

    FILE *fptr;
    fptr = fopen("parallel_universe_coords.txt", "w");
    printf("coords printed! \n");

    for (int i = 0; i <= N; i += J) {
        fprintf(fptr, "%0.2lf \t %0.2lf\n", x[i], y[i]);
    }

    fclose(fptr);
}