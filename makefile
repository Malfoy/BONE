CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time -fmax-errors=3
CFLAGS+=$(LOL)
LDFLAGS=-flto



ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4 -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=Bone

all: $(EXEC)

backbone.o: backbone.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

Bone: backbone.o
	$(CC) -o $@ $^ $(LDFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


