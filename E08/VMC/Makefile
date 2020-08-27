CC = g++
CFLAGS = -Wall -O3 --std=c++11

es8.exe : es8.o random.o
	$(CC) random.o es8.o -o es8.exe
es8.o : es8.cpp es8.h
	$(CC) -c es8.cpp -o es8.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es8.exe
