PARCC = mpicc
CC = gcc
PARRUN = mpirun
CFLAGS = -Ofast -DNDEBUG
M = 1024
N = 1024
F1 = ../out/A.bin
F2 = ../out/B.bin
F3 = ../out/C.bin
NPROC = 27
build: gen calc
run: generate calculate test
generate: gen
	@./gen $(M) $(N) $(F1) $(F2)
gen: gen.o fun.o
	@$(CC) -o gen gen.o fun.o $(CFLAGS)
fun.o: fun.c
	@$(CC) -c -o fun.o fun.c $(CFLAGS)
gen.o: gen.c
	@$(CC) -c -o gen.o gen.c $(CFLAGS)
calc: main.o matrix.o
	@$(PARCC) -o calc main.o matrix.o $(CFLAGS)
main.o: main.c
	@$(PARCC) -c -o main.o main.c $(CFLAGS)
matrix.o: matrix.c
	@$(PARCC) -c -o matrix.o matrix.c $(CFLAGS)
calculate: calc
	@$(PARRUN) -np $(NPROC) ./calc $(F1) $(F2) $(F3)
test:
	@python3 test.py $(F1) $(F2) $(F3)
clean:
	-@rm *.o gen calc $(F1) $(F2) $(F3)
