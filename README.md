# I3
EEIC-2021S-Experiment-I3

# このプログラムについて
学生実験の発展課題のプログラムです。
電力スペクトル密度を利用して、音声の基本周波数を同定します.

特定した基本周波数をピアノの鍵盤上に表示します。

## identify_f0.c
単音の基本周波数を同定します。ウィナーヒンチンの定理を利用しています。

## identify_multi_f0.c
混合音の複数の基本周波数の同定を試みました。完全とは行かないものの、うまく行く場合もあります。

まだ、編集中です。

# 中身
## 鍵盤の表示
```c
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
```
