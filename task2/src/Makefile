CC = gcc
IN = main.c matrix.c
OUT = calc
INCLUDE = -I/usr/local/include
LINK = /usr/local/lib/libpapi.a
CFLAGS = -O0
STDOUT = 
all:
	@date
	@$(CC) $(INCLUDE) $(CFLAGS) $(IN) $(LINK) -o $(OUT)
	@sudo ./$(OUT) $(STDOUT)
	@rm $(OUT)
	@date
