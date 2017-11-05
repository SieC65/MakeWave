FLAGS = -Wall -O3 `root-config --cflags --glibs`

all: MakeWave

MakeWave: main.o MakeWave.o PMT.o
	g++ $(FLAGS) main.o MakeWave.o PMT.o -o MakeWave
	
main.o: main.cpp
	g++ $(FLAGS) -c main.cpp
	
MakeWave.o: MakeWave.cpp
	g++ $(FLAGS) -c MakeWave.cpp
	
PMT.o: PMT.cpp
	g++ $(FLAGS) -c PMT.cpp

clean:
	rm -rf *.o MakeWave
	
#	g++ -c MakeWave.cpp PMT.cpp $(FLAGS) -o MakeWave.o
#	g++ -o MakeWave.exe $(FLAGS) -lrt main.cpp MakeWave.o