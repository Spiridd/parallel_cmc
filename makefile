test:
	@gcc multiplier.c -O3 -o calc
	@for num in 1 2 3 4 5 6; \
	do \
	    ./calc A_bin B_bin C_bin $$num; \
	    echo -n "test $$num: "; \
	    python3 test.py A_bin B_bin C_bin; \
	done
	@rm A_bin B_bin C_bin calc
report:
	@gcc multiplier.c -O3 -o calc
	@./calc A_bin B_bin C_bin 0 times_bin
	@python3 report.py times_bin
	@rm times_bin calc
