CC=	gcc 
CFLAG= -O3 -m64 -w
FLAGS= -lm -lpthread -lrt
HEADER=	util_lib.h fasta.h run_hmm.h
SRCS=	util_lib.c hmm_lib.c run_hmm.c fasta.c 
OBJ=	util_lib.o hmm_lib.o run_hmm.o  fasta.o 
EXEC=  FGS+	

all:  $(OBJ) $(HEADER)
	$(CC)  $(CFLAG) -o $(EXEC) $(OBJ) $(FLAGS)

%.o:%.c
	$(CC) $(CFLAG) -c $< 

clean:
	rm -rf $(OBJ) $(EXEC)

