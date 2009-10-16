OBJS = logtest.o
CC = gcc
CFLAGS = -O3 -lgmp -lmpfr

logtest: $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) -lm

clean:
	rm -f *.o


