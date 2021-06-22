#include <stdio.h>

#define KeyNum 88
void print_keyboard(char *key, int N){
    //int N = 52;
    int i;
    printf("\x1b[30m\x1b[47m|");
    for(i = 0; i < N; i++){
        if(i % 12 == 1 || i % 12 == 4 || i % 12 == 6 || i % 12 == 9 || i % 12 == 11){
            if(key[i] == 1){
                printf("\x1b[41m \x1b[47m");
            }
            else{
                printf("\x1b[40m \x1b[47m");
            }
        }
        else if (i % 12 == 2 || i % 12 == 7){
            printf(" |");
        }
        else{
            printf(" ");
        }
    }
    printf("\n\x1b[47m|");
    for (i = 0; i < N; i++){
        if(i % 12 == 1 || i % 12 == 4 || i % 12 == 6 || i % 12 == 9 || i % 12 == 11){
            printf("|");
        }
        else if (i % 12 == 2 || i % 12 == 7){
            if (key[i] == 1){
                printf("\x1b[41m_\x1b[47m|");            
            }
            else{
                printf("_|");
            }
        }
        else{
            if (key[i] == 1){
                printf("\x1b[41m_\x1b[47m");            
            }
            else{
                printf("_");
            }
        }
    }
    printf("|\x1b[49m\x1b[39m\n");
    fflush(stdout);
}

int main(void){
    char key[KeyNum] = {0};
    key[0] = 1;
    key[1] = 1;
    print_keyboard(key, KeyNum);
    return 0;
}