N_PROC = 10 
N = 100000000
OUT = ./calc
IN = mpi_sieve.c
LINK = -lm
all:
	-@rm out/*.txt
	-@rm $(OUT)
	@mpicc $(IN) $(LINK) -o $(OUT)
	@for num in 1 2 3 4; \
	do \
		echo "run $$num processes"; \
		mpirun -np $$num $(OUT) $(N); \
	done
	@python3 report.py
