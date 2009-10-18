OBJS = gammatest.o
CC = gcc
CFLAGS = -O3 -lgmp -lmpfr

gammatest: $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) -lm

clean:
	rm -f *.o


