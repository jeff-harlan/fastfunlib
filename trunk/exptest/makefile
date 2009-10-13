OBJS = exptest.o
CC = gcc
CFLAGS = -O3 -lgmp -lmpfr

exptest: $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) -lm

clean:
	rm -f *.o


