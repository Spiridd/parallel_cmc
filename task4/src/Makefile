PARCC = mpicc
PARRUN = mpirun
CC = gcc
#CFLAGS = -Ofast
M = 1024
N = 4096
F1 = ../res/A.bin
F2 = ../res/b.bin
F3 = ../res/c.bin
NPROC = 2
# generete - get random matrices
# calucalte - multiply these matrices
# test - check if the result is correct
build: gen calc
run: generate calculate test
generate: gen
	@./gen $(M) $(N) $(F1) $(F2)
gen: gen.o fun.o
	@$(CC) gen.o fun.o -o gen $(CFLAGS)
gen.o: gen.c
	@$(CC) -c -o gen.o gen.c $(CFLAGS)
fun.o: fun.c
	@$(CC) -c -o fun.o fun.c $(CFLAGS)
calculate: calc
	@$(PARRUN) -np $(NPROC) ./calc $(F1) $(F2) $(F3)
calc: main.o matrix.o fun.o
	@$(PARCC) main.o matrix.o fun.o -o calc $(CFLAGS)
main.o: main.c
	@$(PARCC) -c -o main.o main.c $(CFLAGS)
matrix.o: matrix.c
	@$(PARCC) -c -o matrix.o matrix.c $(CFLAGS)
test:
	@python3 test.py $(F1) $(F2) $(F3)
report:
	@python3 report.py
clean:
	-@rm *.o gen calc $(F1) $(F2) $(F3)
