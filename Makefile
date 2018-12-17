CPP = g++
CFLAGS = -Wall -Wextra -Werror -O3 -fPIC -std=c++11 -pthread
INC = -Isrc
LIBS = -lboost_thread -lboost_system

all: memembed mempot

src/pdb.o: src/pdb.cpp
	$(CPP) $(CFLAGS) $(INC) -c src/pdb.cpp -o src/pdb.o

src/ga.o: src/ga.cpp
	$(CPP) $(CFLAGS) $(INC) -c src/ga.cpp -o src/ga.o

src/direct.o: src/direct.cpp
	$(CPP) $(CFLAGS) $(INC) -c src/direct.cpp -o src/direct.o

src/grid.o: src/grid.cpp
	$(CPP) $(CFLAGS) $(INC) -c src/grid.cpp -o src/grid.o

src/main.o: src/main.cpp
	$(CPP) $(CFLAGS) $(INC) -c src/main.cpp -o src/main.o

memembed: src/pdb.o src/ga.o src/grid.o src/direct.o src/main.o
	$(CPP) $(CFLAGS) $(INC) src/main.o src/pdb.o src/ga.o src/grid.o src/direct.o ${LIBS} -o bin/memembed

mempot: src/mempot.cpp
	$(CXX) $(CFLAGS) $(INC) src/mempot.cpp -o bin/mempot

clean:
	rm -f bin/memembed bin/mempot src/*.o examples/*EMBED.pdb

test:
	bin/memembed -t 4,26,46,68 -n out -q 1 -a 2 examples/2x2v.pdb


	
