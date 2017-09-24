test:
	@gcc multiplier.c -O3 -o calc
	@./calc A_bin B_bin C_bin 1
	@echo -n "test 1: "
	@python3 test.py A_bin B_bin C_bin
	@./calc A_bin B_bin C_bin 2 
	@echo -n "test 2: "
	@python3 test.py A_bin B_bin C_bin
	@./calc A_bin B_bin C_bin 3
	@echo -n "test 3: "
	@python3 test.py A_bin B_bin C_bin
	@./calc A_bin B_bin C_bin 4
	@echo -n "test 4: "
	@python3 test.py A_bin B_bin C_bin
	@./calc A_bin B_bin C_bin 5
	@echo -n "test 5: "
	@python3 test.py A_bin B_bin C_bin
	@./calc A_bin B_bin C_bin 6
	@echo -n "test 6: "
	@python3 test.py A_bin B_bin C_bin
	@rm A_bin B_bin C_bin calc
report:
	@gcc multiplier.c -O3 -o calc
	@./calc A_bin B_bin C_bin 0 times_bin
	@python3 report.py times_bin
	@rm times_bin calc
