all:
	g++ -c MakeWave.cpp -Wall -O3 -I. `root-config --cflags --glibs` -o MakeWave.o
	g++ -o MakeWave.exe -Wall -O3 `root-config --cflags --glibs` -lrt main.cpp MakeWave.o