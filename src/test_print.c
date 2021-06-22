#include <stdio.h>

int main(void){
    int N = 52;
    int i;
    printf("\x1b[30m\x1b[47m|");
    for(i = 0; i < N; i++){
        if(i % 7 == 0 || i % 7 == 2 || i % 7 == 3 || i % 7 == 5 || i % 7 == 6 ){
            printf(" \x1b[40m_\x1b[47m");
        }
        else{
            printf(" |");
        }
    }
    printf("\n\x1b[47m");
    for (i = 0; i < N; i++){
        printf("|_");
    }
    printf("|\x1b[49m\x1b[39m \n");
}