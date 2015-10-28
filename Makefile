CC = gcc-5
CFLAGS = -Wall -Wextra
CFLAGS += -O0

PROGS = multic multiz
  
all : $(PROGS)

multic : util.h util.c maf.h maf.c multi_util.h multi_util.c mz_scores.h mz_scores.c mz_preyama.h mz_preyama.c mz_yama.h mz_yama.c charvec.h charvec.c align_util.h align_util.c multic.c
	$(CC) $(CFLAGS) util.c multi_util.c maf.c mz_scores.c charvec.c mz_yama.c mz_preyama.c align_util.c multic.c -o multic

multiz : util.h util.c maf.h maf.c multi_util.h multi_util.c mz_scores.h mz_scores.c mz_preyama.h mz_preyama.c mz_yama.h mz_yama.c charvec.h charvec.c align_util.h align_util.c multiz.c
	$(CC) $(CFLAGS) util.c multi_util.c maf.c mz_scores.c charvec.c mz_yama.c mz_preyama.c align_util.c multiz.c -o multiz

clean:
	$(RM) $(PROGS)
