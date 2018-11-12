FLAGS = -Wall -O3 `root-config --cflags --glibs`

all: MakeWave

MakeWave: main.o MakeWave.o PMT_R11410.o Test.o
	g++ $(FLAGS) main.o MakeWave.o PMT_R11410.o Test.o -o MakeWave

main.o: main.cpp
	g++ $(FLAGS) -c main.cpp

MakeWave.o: MakeWave.cpp
	g++ $(FLAGS) -c MakeWave.cpp

PMT_R11410.o: PMT_R11410.cpp
	g++ $(FLAGS) -c PMT_R11410.cpp

Test.o: Test.cpp
	g++ $(FLAGS) -c Test.cpp

clean:
	rm -rf *.o MakeWave

#	g++ -c MakeWave.cpp PMT.cpp $(FLAGS) -o MakeWave.o
#	g++ -o MakeWave.exe $(FLAGS) -lrt main.cpp MakeWave.o