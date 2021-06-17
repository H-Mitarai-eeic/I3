identify_f0:	./src/identify_f0.c
	gcc -o ./bin/identify_f0 ./src/identify_f0.c -Wall -lm
read_data2:		./src/read_data2.c
	gcc -o ./bin/read_data2 ./src/read_data2.c -Wall -lm