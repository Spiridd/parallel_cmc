CC = gcc
CFLAGS = -O0
INCLUDE = -I/usr/local/include
LINK = /usr/local/lib/libpapi.a
IN = main.c
OUT = run
all:
	@$(CC) $(INCLUDE) $(CFLAGS) $(IN) $(LINK) -o $(OUT)
	@sudo ./$(OUT)
	@rm $(OUT)
