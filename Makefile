####################
# Makefile- for GSL#
# 30-III-09        #
# Nicolas Quesada  #
####################


CC=gcc
CFLAGS=-O4 -Wall -I/usr/include -I.
LFLAGS=-lm -L/usr/lib -lgsl -lgslcblas

%.out:%.o steinhardt.o
	$(CC) $^ $(LFLAGS) -o $@

run:	example.out
	bash runexample.sh

data:
	bash getdata.sh

clean_data:
	rm 13 55 147 4 26 98 38DD 13BB 19AA 75CC 7AA 39AA 39BB 39CC

clean:
	rm -rf *~
	rm -rf *.o
	rm -rf *.out
