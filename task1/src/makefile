CC = gcc
IN = main.c
CFLAGS = -O3
OUT = calc
F1 = A.bin
F2 = B.bin
F3 = C.bin
F4 = times.bin
TYPE = -d
test:
	@$(CC) $(IN) $(CFLAGS) -o $(OUT)
	@for num in 0 1 2 3 4 5; \
	do \
	    ./$(OUT) $(TYPE) $$num $(F1) $(F2) $(F3); \
	    echo -n "test $$num: "; \
	    python3 test.py $(F1) $(F2) $(F3); \
	done
	@rm $(F1) $(F2) $(F3) $(OUT)
report:
	@$(CC) $(IN) $(CFLAGS) -o $(OUT)
	@./$(OUT) $(TYPE) 6 $(F4)
	@python3 report.py $(F4)
	@rm $(F4) $(OUT)
