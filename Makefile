CC=	gcc 
CFLAG= -O3 -m64 -DMAC_SEM
#CFLAG= -O3 -m64
FLAGS= -lm -lpthread  
#FLAGS= -lm -lpthread 
HEADER=	util_lib.h fasta.h
SRCS=	util_lib.c hmm_lib.c run_hmm.c fasta.c
OBJ=	util_lib.o hmm_lib.o run_hmm.o  fasta.o
EXEC=  FGS+	

all:  $(OBJ) $(HEADER)
	$(CC)  $(CFLAG) -o $(EXEC) $(OBJ) $(FLAGS)

%.o:%.c
	$(CC) $(CFLAG) -c $< 

clean:
	rm -rf $(OBJ) $(EXEC)

