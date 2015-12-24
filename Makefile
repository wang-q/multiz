CC = gcc
CFLAGS = -Wall -Wextra
CFLAGS += -O0
LDFLAGS= -L. -lz

PROGS = multiz
  
all : $(PROGS)

multiz : util.h util.c maf.h maf.c multi_util.h multi_util.c mz_scores.h mz_scores.c mz_preyama.h mz_preyama.c mz_yama.h mz_yama.c multiz.c
	$(CC) $(CFLAGS) $(LDFLAGS) util.c multi_util.c maf.c mz_scores.c mz_yama.c mz_preyama.c multiz.c -o multiz

clean:
	$(RM) $(PROGS)
