#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char *argv[]){
    if(argc != 2){
        perror("argc is not 2\n");
        return 1;
    }
    

    char *filename = argv[1]; 
    int fd = open(filename, O_RDONLY, 0644);
    int r_flag, i;
    short data;

    if (fd == -1){
        perror("open");
        exit(1);
    }

    i = 0;
    r_flag = 1;

    while(1){
        r_flag = read(fd, &data, 2);
        if(r_flag !=2){
            break;
        }
        printf("%d %d\n", i, data);
        i++;
    }

    close(fd);
    return 0;
}