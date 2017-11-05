FLAGS = -Wall -O3 `root-config --cflags --glibs`
all:
	g++ -c MakeWave.cpp $(FLAGS) -o MakeWave.o
	g++ -o MakeWave.exe $(FLAGS) -lrt main.cpp MakeWave.o