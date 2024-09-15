#include<stdio.h>
#include <stdlib.h>

double array[6];
#define a 2.1

void output(int length){
    FILE *filepointer = fopen("xyzfile.xyz", "w");
    
    char str[2] = "Si";

    fprintf(filepointer,"%d \n", length);
    fprintf(filepointer,"\n");
    for (int i = 0; i < length ; i++){
        fprintf(filepointer,"%s \t %f \t %f \t %f \n", str, array[i], array[i+1] , array[i+2]);
    }
    fclose(filepointer);
}

int main(){
    //FILE *filepointer = fopen("xyzfile.txt", "w");
    for( int i = 0; i < 6 ; i++ ){
        array[i] = (i+7) + a ;
    }

    int length = 4;
    output(length);
    
  
    return 0;
}