OBJS = logtest2.o
CC = gcc
CFLAGS = -O3 -lgmp -lmpfr

logtest2: $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS) -lm

clean:
	rm -f *.o


