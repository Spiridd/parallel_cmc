PARCC = mpicc
PARRUN = mpirun
CC = gcc
#CFLAGS = -Ofast
M = 1024
N = 4096
F1 = ../out/A.bin
F2 = ../out/b.bin
F3 = ../out/c.bin
NPROC = 2
# generete - get random matrices
# calucalte - multiply these matrices
# test - check if the result is correct
build: gen calc
run: generate calculate test
generate: gen
	@./gen $(M) $(N) $(F1) $(F2)
gen: gen.o fun.o
	@$(PARCC) gen.o fun.o -o gen $(CFLAGS)
fun.o: fun.c
	@$(PARCC) -c -o fun.o fun.c $(CFLAGS)
gen.o: gen.c
	@$(CC) -c -o gen.o gen.c $(CFLAGS)
main.o: main.c
	@$(PARCC) -c -o main.o main.c $(CFLAGS)
calc: main.o fun.o
	@$(PARCC) main.o fun.o -o calc $(CFLAGS)
calculate: calc
	@$(PARRUN) -np $(NPROC) ./calc $(F1) $(F2) $(F3)
test:
	@python3 test.py $(F1) $(F2) $(F3)
clean:
	-@rm *.o gen calc $(F1) $(F2) $(F3)
