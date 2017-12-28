N_PROC = 10 
N = 100000000
OUT = ./calc
IN = mpi_sieve.c
LINK = -lm
all:
	@mpicc $(IN) $(LINK) -o $(OUT)
	@for num in 1 2 3 4; \
	do \
		mpirun -np $$num $(OUT) $(N); \
	done
	@rm $(OUT)
	@python3 report.py
	@rm out/*.txt
